package com.github.servicenow.ds.stats.stl;

/**
 * Encapsulate the complexity of smoothing the cyclic sub-series in a separate class.
 *
 * Created by Jim Crotinger on 13-May-2016.
 */
@SuppressWarnings("WeakerAccess")
public class CyclicSubSeriesSmoother {

	private final double[][] fRawCyclicSubSeries;
	private final double[][] fSmoothedCyclicSubSeries;
	private final double[][] fSubSeriesWeights;

	private final int fPeriodLength;
	private final int fNumPeriods;
	private final int fRemainder;

	private final int fNumPeriodsToExtrapolateBackward;
	private final int fNumPeriodsToExtrapolateForward;

	private final int fWidth;
	private final LoessSmoother.Builder fLoessSmootherFactory;

	/**
	 * Use Builder to simplify complex construction patterns.
	 */
	public static class Builder {
		private int fWidth;
		private int fDegree = 1;
		private int fJump = 1;
		private int fDataLength;
		private int fPeriodicity;
		private int fNumPeriodsBackward = 1;
		private int fNumPeriodsForward = 1;

		public Builder setWidth(int width) {
			fWidth = width;
			return this;
		}

		public Builder setDegree(int degree) {
			if (degree < 0 || degree > 2)
				throw new IllegalArgumentException("Degree must be 0, 1 or 2");

			fDegree = degree;
			return this;
		}

		public Builder setJump(int jump) {
			fJump = jump;
			return this;
		}

		public Builder setDataLength(int dataLength) {
			fDataLength = dataLength;
			return this;
		}

		public Builder setPeriodicity(int periodicity) {
			fPeriodicity = periodicity;
			return this;
		}

		public Builder extrapolateForwardOnly(int periods) {
			fNumPeriodsForward = periods;
			fNumPeriodsBackward = 0;
			return this;
		}

		public Builder setNumPeriodsForward(int periods) {
			fNumPeriodsForward = periods;
			return this;
		}

		public Builder setNumPeriodsBackward(int periods) {
			fNumPeriodsBackward = periods;
			return this;
		}

		public CyclicSubSeriesSmoother build() {
			return new CyclicSubSeriesSmoother(fWidth, fDegree, fJump, fDataLength, fPeriodicity,
					fNumPeriodsBackward, fNumPeriodsForward);
		}
	}

	/**
	 * Create a cyclic sub-series smoother with the specified properties.
	 *
	 * @param width
	 *      - width of the LOESS smoother
	 * @param degree
	 *      - degree of the LOESS smoother
	 * @param jump
	 *      - jump to use in LOESS smoothing
	 * @param dataLength
	 *      - length of the input data
	 * @param periodicity
	 *      - length of the cyclic period
	 * @param numPeriodsToExtrapolateBackward
	 *      - number of periods to extrapolate backward
	 * @param numPeriodsToExtrapolateForward
	 *      - numbers of periods to extrapolate forward
	 */
	CyclicSubSeriesSmoother(int width, int degree, int jump,
							int dataLength, int periodicity,
                            int numPeriodsToExtrapolateBackward, int numPeriodsToExtrapolateForward) {
		fWidth = width;

		fLoessSmootherFactory = new LoessSmoother.Builder().setWidth(width).setJump(jump).setDegree(degree);

		fPeriodLength = periodicity;
		fNumPeriods = dataLength / periodicity;
		fRemainder = dataLength % periodicity;

		fNumPeriodsToExtrapolateBackward = numPeriodsToExtrapolateBackward;
		fNumPeriodsToExtrapolateForward = numPeriodsToExtrapolateForward;

		fRawCyclicSubSeries = new double[periodicity][];
		fSmoothedCyclicSubSeries = new double[periodicity][];
		fSubSeriesWeights = new double[periodicity][];

		// Bookkeeping: Write the data length as
		//
		// n = m * periodicity + r
		//
		// where r < periodicity. The first r sub-series will have length m + 1 and the remaining will have length m.
		// Another way to look at this is that the cycle length is
		//
		// cycleLength = (n - p - 1) / periodicity + 1
		//
		// where p is the index of the cycle that we're currently in.

		for (int period = 0; period < periodicity; ++period) {
			int seriesLength = (period < fRemainder) ? (fNumPeriods + 1) : fNumPeriods;
			fRawCyclicSubSeries[period] = new double[seriesLength];
			fSmoothedCyclicSubSeries[period] = new double[fNumPeriodsToExtrapolateBackward + seriesLength
					+ fNumPeriodsToExtrapolateForward];
			fSubSeriesWeights[period] = new double[seriesLength];
		}
	}

	/**
	 * Run the cyclic sub-series smoother on the specified data, with the specified weights (ignored if null). The
	 * sub-series are reconstructed into a single series in smoothedData.
	 *
	 * @param rawData
	 *            double[] input data
	 * @param smoothedData
	 *            double[] output data
	 * @param weights
	 *            double[] weights to use in the underlying interpolator; ignored if null.
	 */
	public void smoothSeasonal(double[] rawData, double[] smoothedData, double[] weights) {
		extractRawSubSeriesAndWeights(rawData, weights);
		computeSmoothedSubSeries(weights != null);
		reconstructExtendedDataFromSubSeries(smoothedData);
		// SeasonalTrendLoess.dumpDebugData("extended seasonal", smoothedData);
	}

	private void computeSmoothedSubSeries(boolean useResidualWeights) {
		for (int period = 0; period < fPeriodLength; ++period) {
			double[] weights = useResidualWeights ? fSubSeriesWeights[period] : null;
			double[] rawData = fRawCyclicSubSeries[period];
			double[] smoothedData = fSmoothedCyclicSubSeries[period];

			smoothOneSubSeries(weights, rawData, smoothedData);

			// dumpCyclicSubseriesDebugData(period, rawData.length, smoothedData, rawData);
		}
	}

	private void extractRawSubSeriesAndWeights(double[] data, double[] weights) {
		for (int period = 0; period < fPeriodLength; ++period) {
			final int cycleLength = (period < fRemainder) ? (fNumPeriods + 1) : fNumPeriods;
			for (int i = 0; i < cycleLength; ++i) {
				fRawCyclicSubSeries[period][i] = data[i * fPeriodLength + period];
				if (weights != null) {
					fSubSeriesWeights[period][i] = weights[i * fPeriodLength + period];
				}
			}
		}
	}

	private void reconstructExtendedDataFromSubSeries(double[] data) {
		// Copy this smoothed cyclic sub-series to the extendedSeasonal work array.
		for (int period = 0; period < fPeriodLength; ++period) {
			final int cycleLength = (period < fRemainder) ? (fNumPeriods + 1) : fNumPeriods;
			for (int i = 0; i < fNumPeriodsToExtrapolateBackward + cycleLength + fNumPeriodsToExtrapolateForward; ++i) {
				data[i * fPeriodLength + period] = fSmoothedCyclicSubSeries[period][i];
			}
		}
	}

	/**
	 * Use LOESS interpolation on each of the cyclic sub-series (e.g. in monthly data, smooth the Januaries, Februaries,
	 * etc.).
	 *
	 * @param weights
	 * @param rawData
	 * @param smoothedData
	 */
	 private void smoothOneSubSeries(double[] weights, double[] rawData, double[] smoothedData) {

		final int cycleLength = rawData.length;

		// Smooth the cyclic sub-series with LOESS and then extrapolate one place beyond each end.

		LoessSmoother smoother = fLoessSmootherFactory.setData(rawData).setExternalWeights(weights).build();

		// Copy, shifting by 1 to leave room for the extrapolated point at the beginning.

		System.arraycopy(smoother.smooth(), 0, smoothedData, fNumPeriodsToExtrapolateBackward, cycleLength);

		LoessInterpolator interpolator = smoother.getInterpolator();

		// Extrapolate from the leftmost "width" points to the "-1" position
		int left = 0;
		int right = left + fWidth - 1;
		right = Math.min(right, cycleLength - 1);
		int leftValue = fNumPeriodsToExtrapolateBackward;

		for (int i = 1; i <= fNumPeriodsToExtrapolateBackward; i++) {
			Double ys = interpolator.smoothOnePoint(-i, left, right);
			smoothedData[leftValue - i] = ys == null ? smoothedData[leftValue] : ys;
		}

		// Extrapolate from the rightmost "width" points to the "length" position (one past the array end).
		right = cycleLength - 1;
		left = right - fWidth + 1;
		left = Math.max(0, left);
		int rightValue = fNumPeriodsToExtrapolateBackward + right;

		for (int i = 1; i <= fNumPeriodsToExtrapolateForward; i++) {
			Double ys = interpolator.smoothOnePoint(right + i, left, right);
			smoothedData[rightValue + i] = ys == null ? smoothedData[rightValue] : ys;
		}
	}

//	@SuppressWarnings("unused")
//	private static void dumpCyclicSubseriesDebugData(int p, int cycleLength, double[] smoothedData, double[] inputData) {
//		System.out.println(String.format("subcycle smoother at j = %d", p));
//		System.out.println(String.format("                              , smoothed(%2d) = %22.15e", 0, smoothedData[0]));
//		for (int i = 0; i < cycleLength; ++i) {
//			System.out.println(String.format("y(%2d) = %22.15e, smoothed(%2d) = %22.15e", i, inputData[i], i + 1,
//			smoothedData[i + 1]));
//		}
//		System.out.println(String.format("                              , smoothed(%2d) = %22.15e", cycleLength + 1,
//		smoothedData[cycleLength + 1]));
//	}
}

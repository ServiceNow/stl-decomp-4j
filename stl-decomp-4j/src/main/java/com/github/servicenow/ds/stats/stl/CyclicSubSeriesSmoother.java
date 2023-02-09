package com.github.servicenow.ds.stats.stl;

import java.util.Arrays;

/**
 * Encapsulate the complexity of smoothing the cyclic sub-series in a separate class.
 * <p>
 * Created by Jim Crotinger on 13-May-2016.
 */
@SuppressWarnings("WeakerAccess")
public class CyclicSubSeriesSmoother {

	private final double[][] fRawCyclicSubSeries;
	private final double[][] fSmoothedCyclicSubSeries;

	// reshaped version of exogenous inputs with an extra dimension created to form subseries with respect to the periodicity,
	// thereby leading to the dimensions [periodicity][numOfExogData][numOfDataInEachCycle]
	private final double[][][] fExogenousCyclicSeries;

	// controls the output of non-exogenous component from the fit
	private final boolean fOutputNonExogenousPart;

	private final double[][] fSubSeriesWeights;

	private final int fPeriodLength;
	private final int fNumPeriods;
	private final int fRemainder;

	private final int fNumPeriodsToExtrapolateBackward;
	private final int fNumPeriodsToExtrapolateForward;

	private final int fWidth;
	private final int fDegree;
	private final LoessSmoother.Builder fLoessSmootherFactory;

	/**
	 * Use Builder to simplify complex construction patterns.
	 */
	public static class Builder {
		private Integer fWidth = null;
		private Integer fDataLength = null;
		private Integer fPeriodicity = null;
		private Integer fNumPeriodsBackward = null;
		private Integer fNumPeriodsForward = null;
		private int fDegree = 1;
		private int fJump = 1;
		private int fNumExogenousInputs;
		private boolean fOutputNonExogenousPart;

		/**
		 * Set the width of the LOESS smoother used to smooth each seasonal sub-series.
		 *
		 * @param width width of the LOESS smoother
		 * @return this
		 */
		public Builder setWidth(int width) {
			fWidth = width;
			return this;
		}

		/**
		 * Set the degree of the LOESS smoother used to smooth each seasonal sub-series.
		 *
		 * @param degree degree of the LOESS smoother
		 * @return this
		 */
		public Builder setDegree(int degree) {
			if (degree < 0 || degree > 2)
				throw new IllegalArgumentException("Degree must be 0, 1 or 2");

			fDegree = degree;
			return this;
		}

		/**
		 * Set the jump (number of points to skip) between LOESS interpolations when smoothing the seasonal sub-series.
		 * <p>
		 * Defaults to 1 (computes LOESS interpolation at each point).
		 *
		 * @param jump jump (number of points to skip) in the LOESS smoother
		 * @return this
		 */
		public Builder setJump(int jump) {
			fJump = jump;
			return this;
		}

		/**
		 * Set the total length of the data that will be deconstructed into cyclic sub-series.
		 *
		 * @param dataLength total length of the data
		 * @return this
		 */
		public Builder setDataLength(int dataLength) {
			fDataLength = dataLength;
			return this;
		}

		/**
		 * Set the period of the data's seasonality.
		 *
		 * @param periodicity number of data points in each season or period
		 * @return this
		 */
		public Builder setPeriodicity(int periodicity) {
			fPeriodicity = periodicity;
			return this;
		}

		/**
		 * Set the number of exogenous inputs
		 *
		 * @param numExogenousInputs number of exogenous inputs
		 * @return this
		 */
		public Builder setNumExogenousInputs(int numExogenousInputs) {
			fNumExogenousInputs = numExogenousInputs;
			return this;
		}

		/**
		 * Set the boolean for outputting only the nonexogenous part from the smoother.
		 *
		 * @param outputNonExogenousPart boolean if True then only outputs the const+trend from the smoother.
		 * @return this
		 */
		public Builder setOutputNonExogenousPart(boolean outputNonExogenousPart) {
			fOutputNonExogenousPart = outputNonExogenousPart;
			return this;
		}

		/**
		 * Construct a smoother that will extrapolate forward only by the specified number of periods.
		 *
		 * @param periods number of periods to extrapolate
		 * @return this
		 */
		public Builder extrapolateForwardOnly(int periods) {
			fNumPeriodsForward = periods;
			fNumPeriodsBackward = 0;
			return this;
		}

		/**
		 * Construct a smoother that extrapolates forward and backward by the specified number of periods.
		 *
		 * @param periods number of periods to extrapolate
		 * @return this
		 */
		public Builder extrapolateForwardAndBack(int periods) {
			fNumPeriodsForward = periods;
			fNumPeriodsBackward = periods;
			return this;
		}

		/**
		 * Set the number of periods to extrapolate forward.
		 * <p>
		 * Defaults to 1.
		 *
		 * @param periods number of periods to extrapolate
		 * @return this
		 */
		public Builder setNumPeriodsForward(int periods) {
			fNumPeriodsForward = periods;
			return this;
		}

		/**
		 * Set the number of periods to extrapolate backward.
		 * <p>
		 * Defaults to 1.
		 *
		 * @param periods number of periods to extrapolate
		 * @return this
		 */
		public Builder setNumPeriodsBackward(int periods) {
			fNumPeriodsBackward = periods;
			return this;
		}

		/**
		 * Build the sub-series smoother.
		 *
		 * @return new CyclicSubSeriesSmoother
		 */
		public CyclicSubSeriesSmoother build() {
			checkSanity();

			return new CyclicSubSeriesSmoother(fWidth, fDegree, fJump, fDataLength, fPeriodicity,
					fNumPeriodsBackward, fNumPeriodsForward, fNumExogenousInputs, fOutputNonExogenousPart);
		}

		private void checkSanity() {
			if (fWidth == null)
				throw new IllegalArgumentException(
						"CyclicSubSeriesSmoother.Builder: setWidth must be called before building the smoother.");

			if (fPeriodicity == null)
				throw new IllegalArgumentException(
						"CyclicSubSeriesSmoother.Builder: setPeriodicity must be called before building the smoother.");

			if (fDataLength == null)
				throw new IllegalArgumentException(
						"CyclicSubSeriesSmoother.Builder: setDataLength must be called before building the smoother.");

			if (fNumPeriodsBackward == null || fNumPeriodsForward == null)
				throw new IllegalArgumentException(
						"CyclicSubSeriesSmoother.Builder: Extrapolation settings must be provided.");
		}
	}

	/**
	 * Create a cyclic sub-series smoother with the specified properties.
	 *
	 * @param width                           width of the LOESS smoother
	 * @param degree                          degree of the LOESS smoother
	 * @param jump                            jump to use in LOESS smoothing
	 * @param dataLength                      length of the input data
	 * @param periodicity                     length of the cyclic period
	 * @param numPeriodsToExtrapolateBackward number of periods to extrapolate backward
	 * @param numPeriodsToExtrapolateForward  numbers of periods to extrapolate forward
	 * @param numExogenousInputs              number of exogenous inputs
	 * @param outputNonExogenousPart               if true it separates trend from exogenous inputs and outputs the former
	 *                                        if false it outputs the whole (smoothed) fit to the trend and exogenous inputs
	 */
	CyclicSubSeriesSmoother(int width, int degree, int jump,
	                        int dataLength, int periodicity,
	                        int numPeriodsToExtrapolateBackward, int numPeriodsToExtrapolateForward,
	                        int numExogenousInputs, boolean outputNonExogenousPart) {
		fWidth = width;
		fDegree = degree;

		fLoessSmootherFactory = new LoessSmoother.Builder().setWidth(width).setJump(jump).setDegree(degree);

		// For overdetermined system in the smoothing one should have (dataLength/periodicity) >= (numExogenousInputs + degree + 1)
		// This can be made more strict via multiplying the rhs. by 2-5 or adding some integer where the latter is done below.
		if (outputNonExogenousPart && ((dataLength/periodicity) < (numExogenousInputs + degree + 4)) )
			periodicity = 1; // if periodicity > 1 does not allow for subseries fit then fit the whole series at once.

		fPeriodLength = periodicity;
		fNumPeriods = dataLength / periodicity;
		fRemainder = dataLength % periodicity;

		fOutputNonExogenousPart = outputNonExogenousPart;
		fNumPeriodsToExtrapolateBackward = numPeriodsToExtrapolateBackward;
		fNumPeriodsToExtrapolateForward = numPeriodsToExtrapolateForward;

		fRawCyclicSubSeries = new double[periodicity][];
		fSmoothedCyclicSubSeries = new double[periodicity][];
		fExogenousCyclicSeries = new double[periodicity][][];
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
			fExogenousCyclicSeries[period] = null;
			fSmoothedCyclicSubSeries[period] = new double[fNumPeriodsToExtrapolateBackward + seriesLength
					+ fNumPeriodsToExtrapolateForward];
			fSubSeriesWeights[period] = new double[seriesLength];
		}
	}

	/**
	 * Run the cyclic sub-series smoother on the specified data, with the specified weights (ignored if null). The
	 * sub-series are reconstructed into a single series in smoothedData.
	 *
	 * @param rawData      input data
	 * @param smoothedData output data
	 * @param exogenousinputs exogenous inputs data where each row is an input
	 * @param weights      weights to use in the underlying interpolator; ignored if null.
	 */
	public void smoothSeasonal(double[] rawData, double[] smoothedData, double[][] exogenousinputs, double[] weights) {
		if (exogenousinputs == null)
			extractRawSubSeriesAndWeights(rawData, weights);
		else
			extractRawSubSeriesAndWeights(rawData, exogenousinputs, weights);
		computeSmoothedSubSeries(weights != null);
		reconstructExtendedDataFromSubSeries(smoothedData);
		// SeasonalTrendLoess.dumpDebugData("extended seasonal", smoothedData);
	}

	public void smoothSeasonal(double[] rawData, double[] smoothedData, double[] weights) {
		smoothSeasonal(rawData, smoothedData, null, weights);
	}

	private void computeSmoothedSubSeries(boolean useResidualWeights) {
		for (int period = 0; period < fPeriodLength; ++period) {
			double[] weights = useResidualWeights ? fSubSeriesWeights[period] : null;
			double[] rawData = fRawCyclicSubSeries[period];
			double[] smoothedData = fSmoothedCyclicSubSeries[period];
			double[][] exogenousCyclicSeries = fExogenousCyclicSeries[period];

			smoothOneSubSeries(weights, rawData, smoothedData, exogenousCyclicSeries);

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

	private void extractRawSubSeriesAndWeights(double[] data, double[][] exogenousinputs, double[] weights) {
		for (int period = 0; period < fPeriodLength; ++period) {
			fExogenousCyclicSeries[period] = new double[exogenousinputs.length][fRawCyclicSubSeries[period].length];
			final int cycleLength = (period < fRemainder) ? (fNumPeriods + 1) : fNumPeriods;
			for (int i = 0; i < cycleLength; ++i) {
				int index = i * fPeriodLength + period;
				fRawCyclicSubSeries[period][i] = data[index];
				for (int j = 0; j < exogenousinputs.length; ++j) {
					fExogenousCyclicSeries[period][j][i] = exogenousinputs[j][index];
				}
				if (weights != null) {
					fSubSeriesWeights[period][i] = weights[index];
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
	 * @param weights      external weights for interpolation
	 * @param rawData      input data to be smoothed
	 * @param smoothedData output smoothed data
	 */
	private void smoothOneSubSeries(double[] weights, double[] rawData, double[] smoothedData, double[][] exogenousCyclicSeries) {

		final int cycleLength = rawData.length;

		// Smooth the cyclic sub-series with LOESS and then extrapolate one place beyond each end.

		LoessSmoother smoother = fLoessSmootherFactory
				.setData(rawData)
				.setExogenousInputs(exogenousCyclicSeries)
				.setOutputNonExogenousPart(fOutputNonExogenousPart)
				.setExternalWeights(weights)
				.build();
				
		// Copy, shifting by 1 to leave room for the extrapolated point at the beginning.

		System.arraycopy(smoother.smooth(), 0, smoothedData, fNumPeriodsToExtrapolateBackward, cycleLength);

		LoessInterpolator interpolator;
		if (fNumPeriodsToExtrapolateForward > 0 && fOutputNonExogenousPart)
			interpolator = new LoessInterpolator.Builder()
					.setWidth(fWidth)
					.setDegree(fDegree)
					.setOutputNonExogenousPart(fOutputNonExogenousPart)
					.setExternalWeights(weights)
					.interpolate(Arrays.copyOf(smoothedData, rawData.length), null); // if with exogs, remove them and forecast only non-exog part.
		else
			interpolator = smoother.getInterpolator();

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

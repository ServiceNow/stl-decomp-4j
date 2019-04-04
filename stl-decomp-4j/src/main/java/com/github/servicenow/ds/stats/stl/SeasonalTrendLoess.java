package com.github.servicenow.ds.stats.stl;

import java.util.Arrays;

import static com.github.servicenow.ds.stats.TimeSeriesUtilities.simpleMovingAverage;

/**
 * Java implementation of the Seasonal-Trend-Loess algorithm for evenly spaced data. This is basically a direct port of
 * the RATFOR from the netlib stl package.
 * <p>
 * Created by Jim Crotinger on 18-Apr-2016.
 */
public class SeasonalTrendLoess {

	private final double[] fData;

	private Decomposition fDecomposition;

	private final int fPeriodLength;
	private final LoessSettings fSeasonalSettings;
	private final LoessSettings fTrendSettings;
	private final LoessSettings fLowpassSettings;
	private final int fInnerIterations;
	private final int fRobustIterations;

	private final double[] fDetrend;
	private final double[] fExtendedSeasonal;

	private double[] fDeSeasonalized; // TODO: Garbage - can this be made in-place?

	private final CyclicSubSeriesSmoother fCyclicSubSeriesSmoother;
	private final LoessSmoother.Builder fLoessSmootherFactory;
	private final LoessSmoother.Builder fLowpassLoessFactory;

	/**
	 * Builder class for SeasonalTrendLoess decomposition
	 */
	public static class Builder {
		private Integer fPeriodLength = null;

		private Integer fSeasonalWidth = null;
		private Integer fSeasonalJump = null;
		private Integer fSeasonalDegree = null;

		private Integer fTrendWidth = null;
		private Integer fTrendJump = null;
		private Integer fTrendDegree = null;

		private Integer fLowpassWidth = null;
		private Integer fLowpassJump = null;
		private int fLowpassDegree = 1;

		// Following the R interface, we default to "non-robust"

		private int fInnerIterations = 2;
		private int fRobustIterations = 0;

		// Following the R interface, we implement a "periodic" flag that defaults to false.

		private boolean fPeriodic = false;
		private boolean fFlatTrend = false;
		private boolean fLinearTrend = false;

		private LoessSettings buildSettings(int width, int degree, Integer jump) {
			if (jump == null) {
				return new LoessSettings(width, degree);
			} else {
				return new LoessSettings(width, degree, jump);
			}
		}

		/**
		 * Set the period length for the STL seasonal decomposition.
		 * <p>
		 * Required - no default.
		 *
		 * @param period period length (number of data points in each season or period)
		 * @return this
		 */
		public Builder setPeriodLength(int period) {
			if (period < 2)
				throw new IllegalArgumentException("periodicity must be at least 2");

			fPeriodLength = period;
			return this;
		}

		/**
		 * Set the LOESS width (in data points) used to smooth the seasonal sub-series.
		 * <p>
		 * Required unless setPeriodic is called.
		 *
		 * @param width LOESS width for the seasonal sub-series
		 * @return this
		 */
		public Builder setSeasonalWidth(int width) {
			fSeasonalWidth = width;
			return this;
		}

		/**
		 * Set the LOESS degree used to smooth the seasonal sub-series.
		 * <p>
		 * Defaults to 1.
		 *
		 * @param degree LOESS degree for the seasonal sub-series
		 * @return this
		 */
		public Builder setSeasonalDegree(int degree) {
			fSeasonalDegree = degree;
			return this;
		}

		/**
		 * Set the jump (number of points to skip) between LOESS interpolations when smoothing the seasonal sub-series.
		 * <p>
		 * Defaults to 10% of the smoother width.
		 *
		 * @param jump LOESS jump (number of points to skip) for the seasonal sub-series
		 * @return this
		 */
		public Builder setSeasonalJump(int jump) {
			fSeasonalJump = jump;
			return this;
		}

		/**
		 * Set the LOESS width (in data points) used to smooth the trend.
		 * <p>
		 * Defaults to (int) (1.5 * periodLength / (1 - 1.5 / seasonalWidth) + 0.5)
		 *
		 * @param width LOESS with for the trend component
		 * @return this
		 */
		public Builder setTrendWidth(int width) {
			fTrendWidth = width;
			return this;
		}

		/**
		 * Set the LOESS degree used to smooth the trend.
		 * <p>
		 * Defaults to 1.
		 *
		 * @param degree LOESS degree for the trend component
		 * @return this
		 */
		public Builder setTrendDegree(int degree) {
			fTrendDegree = degree;
			return this;
		}

		/**
		 * Set the jump (number of points to skip) between LOESS interpolations used when smoothing the trend.
		 * <p>
		 * Defaults to 10% of the smoother width.
		 *
		 * @param jump LOESS jump (number of points to skip) for the trend component
		 * @return this
		 */
		public Builder setTrendJump(int jump) {
			fTrendJump = jump;
			return this;
		}

		/**
		 * Set the LOESS width (in data points) used by the low-pass filter step.
		 * <p>
		 * Defaults to the period length.
		 *
		 * @param width LOESS width for the low-pass step
		 * @return this
		 */
		public Builder setLowpassWidth(int width) {
			fLowpassWidth = width;
			return this;
		}

		/**
		 * Set the LOESS degree used by the low-pass filter step.
		 * <p>
		 * Defaults to 1.
		 *
		 * @param degree LOESS degree for the low-pass step
		 * @return this
		 */
		public Builder setLowpassDegree(int degree) {
			fLowpassDegree = degree;
			return this;
		}

		/**
		 * Set the jump (number of points to skip) between LOESS interpolations used by the low-pass filter step.
		 * <p>
		 * Defaults to 10% of the smoother width.
		 *
		 * @param jump LOESS jump (number of points to skip) for the low-pass step
		 * @return this
		 */
		public Builder setLowpassJump(int jump) {
			fLowpassJump = jump;
			return this;
		}

		/**
		 * Set the number of STL inner iterations.
		 * <p>
		 * Required, but also set by setRobust, setNonRobust, setRobustFlag.
		 *
		 * @param ni number of inner iterations
		 * @return this
		 */
		public Builder setInnerIterations(int ni) {
			fInnerIterations = ni;
			return this;
		}

		/**
		 * Set the number of STL robustness (outer) iterations.
		 * <p>
		 * Required, but also set by setRobust, setNonRobust, setRobustFlag.
		 *
		 * @param no number of outer iterations
		 * @return this
		 */
		public Builder setRobustnessIterations(int no) {
			fRobustIterations = no;
			return this;
		}

		/**
		 * Set the default robust STL iteration counts (15 robustness iterations, 1 inner iteration).
		 *
		 * @return this
		 */
		public Builder setRobust() {
			fInnerIterations = 1;
			fRobustIterations = 15;
			return this;
		}

		/**
		 * Set the default non-robust STL iteration counts (0 robustness iterations, 2 inner iterations).
		 *
		 * @return this
		 */
		public Builder setNonRobust() {
			fInnerIterations = 2;
			fRobustIterations = 0;
			return this;
		}

		/**
		 * Set the robustness according to a flag; e.g. setRobust if true, setNonRobust if false.
		 *
		 * @param robust true to be robust
		 * @return this
		 */
		public Builder setRobustFlag(boolean robust) {
			return robust ? setRobust() : setNonRobust();
		}

		/**
		 * Constrain the seasonal component to be exactly periodic.
		 *
		 * @return this
		 */
		public Builder setPeriodic() {
			fPeriodic = true;
			return this;
		}

		/**
		 * Set the trend smoother force a flat trend. (Degree == 0, Large Loess Width)
		 *
		 * @return this
		 */
		public Builder setFlatTrend() {
			fLinearTrend = false;
			fFlatTrend = true;
			return this;
		}

		/**
		 * Set the trend smoother force a linear trend. (Degree == 1, Large Loess Width)
		 *
		 * @return this
		 */
		public Builder setLinearTrend() {
			fFlatTrend = false;
			fLinearTrend = true;
			return this;
		}

		/**
		 * Construct the smoother.
		 *
		 * @param data the data to be smoothed
		 * @return a new SeasonalTrendLoess object
		 */
		public SeasonalTrendLoess buildSmoother(double[] data) {
			sanityCheck(data);

			if (fPeriodic) {
				fSeasonalWidth = 100 * data.length;
				fSeasonalDegree = 0;
			} else if (fSeasonalDegree == null) {
				fSeasonalDegree = 1;
			}

			LoessSettings seasonalSettings = buildSettings(fSeasonalWidth, fSeasonalDegree, fSeasonalJump);

			if (fFlatTrend) {
				fTrendWidth = 100 * fPeriodLength * data.length;
				fTrendDegree = 0;
			} else if (fLinearTrend) {
				fTrendWidth = 100 * fPeriodLength * data.length;
				fTrendDegree = 1;
			} else if (fTrendDegree == null) {
				fTrendDegree = 1;
			}

			if (fTrendWidth == null)
				fTrendWidth = calcDefaultTrendWidth(fPeriodLength, fSeasonalWidth);

			LoessSettings trendSettings = buildSettings(fTrendWidth, fTrendDegree, fTrendJump);

			if (fLowpassWidth == null)
				fLowpassWidth = fPeriodLength;

			LoessSettings lowpassSettings = buildSettings(fLowpassWidth, fLowpassDegree, fLowpassJump);

			return new SeasonalTrendLoess(data, fPeriodLength, fInnerIterations, fRobustIterations,
					seasonalSettings, trendSettings, lowpassSettings);
		}

		private static int calcDefaultTrendWidth(int periodicity, int seasonalWidth) {
			// This formula is based on a numerical stability analysis in the original paper.
			return (int) (1.5 * periodicity / (1 - 1.5 / seasonalWidth) + 0.5);
		}

		private void sanityCheck(double[] data) {
			if (data == null)
				throw new IllegalArgumentException(
						"SeasonalTrendLoess.Builder: Data array must be non-null");

			if (fPeriodLength == null)
				throw new IllegalArgumentException(
						"SeasonalTrendLoess.Builder: Period Length must be specified");

			if (data.length < 2 * fPeriodLength)
				throw new IllegalArgumentException(
						"SeasonalTrendLoess.Builder: Data series must be at least 2 * periodicity in length");

			if (fPeriodic) {
				int massiveWidth = 100 * data.length;

				boolean periodicConsistent =
						fSeasonalDegree != null && fSeasonalWidth != null &&
						fSeasonalWidth ==  massiveWidth && fSeasonalDegree == 0;

				if (fSeasonalWidth != null && !periodicConsistent)
					throw new IllegalArgumentException(
							"SeasonalTrendLoess.Builder: setSeasonalWidth and setPeriodic cannot both be called.");

				if (fSeasonalDegree != null && !periodicConsistent)
					throw new IllegalArgumentException(
							"SeasonalTrendLoess.Builder: setSeasonalDegree and setPeriodic cannot both be called.");

				if (fSeasonalJump != null)
					throw new IllegalArgumentException(
							"SeasonalTrendLoess.Builder: setSeasonalJump and setPeriodic cannot both be called.");
			} else {
				if (fSeasonalWidth == null)
					throw new IllegalArgumentException(
							"SeasonalTrendLoess.Builder: setSeasonalWidth or setPeriodic must be called.");
			}

			if (fFlatTrend) {

				int massiveWidth = 100 * fPeriodLength * data.length;

				boolean flatTrendConsistent = fTrendWidth != null && fTrendDegree != null &&
						fTrendWidth == massiveWidth && fTrendDegree == 0;

				if (fTrendWidth != null && !flatTrendConsistent)
					throw new IllegalArgumentException(
							"SeasonalTrendLoess.Builder: setTrendWidth incompatible with flat trend.");

				if (fTrendDegree != null && !flatTrendConsistent)
					throw new IllegalArgumentException(
							"SeasonalTrendLoess.Builder: setTrendDegree incompatible with flat trend.");

				if (fTrendJump != null)
					throw new IllegalArgumentException(
							"SeasonalTrendLoess.Builder: setTrendJump incompatible with flat trend.");
			}

			if (fLinearTrend) {

				int massiveWidth = 100 * fPeriodLength * data.length;

				boolean linearTrendConsistent = fTrendWidth != null && fTrendDegree != null &&
						fTrendWidth == massiveWidth && fTrendDegree == 1;

				if (fTrendWidth != null && !linearTrendConsistent)
					throw new IllegalArgumentException(
							"SeasonalTrendLoess.Builder: setTrendWidth incompatible with linear trend.");

				if (fTrendDegree != null && !linearTrendConsistent)
					throw new IllegalArgumentException(
							"SeasonalTrendLoess.Builder: setTrendDegree incompatible with linear trend.");

				if (fTrendJump != null)
					throw new IllegalArgumentException(
							"SeasonalTrendLoess.Builder: setTrendJump incompatible with linear trend.");
			}
		}
	}

	/**
	 * Construct STL specifying full details of the LOESS smoothers via LoessSettings objects.
	 *
	 * @param data             the data to be decomposed
	 * @param periodicity      the periodicity of the data
	 * @param ni               the number of inner iterations
	 * @param no               the number of outer "robustness" iterations
	 * @param seasonalSettings the settings for the LOESS smoother for the cyclic sub-series
	 * @param trendSettings    the settings for the LOESS smoother for the trend component
	 * @param lowpassSettings  the settings for the LOESS smoother used in de-seasonalizing
	 */
	// Could be private but causes a hidden class to be generated in order for the Builder to have access.
	@SuppressWarnings("WeakerAccess")
	SeasonalTrendLoess(double[] data, int periodicity, int ni, int no, LoessSettings seasonalSettings,
	                   LoessSettings trendSettings, LoessSettings lowpassSettings) {

		fData = data;

		final int size = data.length;

		fPeriodLength = periodicity;
		fSeasonalSettings = seasonalSettings;
		fTrendSettings = trendSettings;
		fLowpassSettings = lowpassSettings;
		fInnerIterations = ni;
		fRobustIterations = no;

		fLoessSmootherFactory = new LoessSmoother.Builder() //
				.setWidth(fTrendSettings.getWidth()) //
				.setDegree(fTrendSettings.getDegree()) //
				.setJump(fTrendSettings.getJump());

		fLowpassLoessFactory = new LoessSmoother.Builder() //
				.setWidth(fLowpassSettings.getWidth()) //
				.setDegree(fLowpassSettings.getDegree()) //
				.setJump(fLowpassSettings.getJump());

		fCyclicSubSeriesSmoother = new CyclicSubSeriesSmoother.Builder() //
				.setWidth(seasonalSettings.getWidth()) //
				.setDegree(seasonalSettings.getDegree()) //
				.setJump(seasonalSettings.getJump()) //
				.setDataLength(size) //
				.extrapolateForwardAndBack(1) //
				.setPeriodicity(periodicity).build();

		fDetrend = new double[size];
		fExtendedSeasonal = new double[size + 2 * fPeriodLength];
	}

	/**
	 * Factory method to perform a non-robust STL decomposition enforcing strict periodicity.
	 * <p>
	 * Meant for diagnostic purposes only.
	 *
	 * @param data        the data to analyze
	 * @param periodicity the (suspected) periodicity of the data
	 * @return SeasonalTrendLoess object with the decomposition already performed.
	 */
	public static Decomposition performPeriodicDecomposition(
			double[] data,
			int periodicity
	) {
		// The LOESS interpolator with degree 0 and a very long window (arbitrarily chosen to be 100 times the length of
		// the array) will interpolate all points as the average value of the series. This particular setting is used
		// for smoothing the seasonal sub-cycles, so the end result is that the seasonal component of the decomposition
		// is exactly periodic.

		// This fit is for diagnostic purposes, so we just do a single inner iteration.

		SeasonalTrendLoess stl = new SeasonalTrendLoess.Builder()
				.setPeriodLength(periodicity) //
				.setSeasonalWidth(100 * data.length) //
				.setSeasonalDegree(0) //
				.setInnerIterations(1) //
				.setRobustnessIterations(0) //
				.buildSmoother(data);

		return stl.decompose();
	}

	/**
	 * Factory method to perform a (somewhat) robust STL decomposition enforcing strict periodicity.
	 * <p>
	 * Meant for diagnostic purposes only.
	 *
	 * @param data        the data to analyze
	 * @param periodicity the (suspected) periodicity of the data
	 * @return SeasonalTrendLoess object with the decomposition already performed.
	 */
	public static Decomposition performRobustPeriodicDecomposition(
			double[] data,
			int periodicity
	) {
		// The LOESS interpolator with degree 0 and a very long window (arbitrarily chosen to be 100 times the length of
		// the array) will interpolate all points as the average value of the series. This particular setting is used
		// for smoothing the seasonal sub-cycles, so the end result is that the seasonal component of the decomposition
		// is exactly periodic.

		// This fit is for diagnostic purposes, so we just do a single inner and outer iteration.

		SeasonalTrendLoess stl = new SeasonalTrendLoess.Builder()
				.setPeriodLength(periodicity) //
				.setSeasonalWidth(100 * data.length) //
				.setSeasonalDegree(0) //
				.setInnerIterations(1) //
				.setRobustnessIterations(1) //
				.buildSmoother(data);

		return stl.decompose();
	}

	/**
	 * Simple class to hold the results of the STL decomposition.
	 */
	public static class Decomposition {

		private final double[] fData;
		private final double[] fTrend;
		private final double[] fSeasonal;
		private final double[] fResiduals;
		private final double[] fWeights;

		/**
		 * Initialize Decomposition object from the original data.
		 * <p>
		 * Allocates space for the decomposition and initializes the weights to 1.
		 *
		 * @param data input data
		 */
		Decomposition(double[] data) {
			fData = data;

			int size = fData.length;
			fTrend = new double[size];
			fSeasonal = new double[size];
			fResiduals = new double[size];
			fWeights = new double[size];

			Arrays.fill(fWeights, 1.0);
		}

		/**
		 * Get the data array that was decomposed
		 *
		 * @return double[] the original data
		 */
		public double[] getData() {
			return fData;
		}

		/**
		 * Get the trend component of the decomposition
		 *
		 * @return double[] the trend component
		 */
		public double[] getTrend() {
			return fTrend;
		}

		/**
		 * Get the seasonal component of the decomposition
		 *
		 * @return double[] the seasonal component
		 */
		public double[] getSeasonal() {
			return fSeasonal;
		}

		/**
		 * Get the residual remaining after removing the seasonality and trend from the data.
		 *
		 * @return double[] the residuals
		 */
		public double[] getResidual() {
			return fResiduals;
		}

		/**
		 * Get the robustness weights used in the calculation. Places where the weights are near zero indicate outliers that
		 * were effectively ignored during the decomposition. (Only applicable if robustness iterations are performed.)
		 *
		 * @return double[] the robustness weights
		 */
		public double[] getWeights() {
			return fWeights;
		}

		private void updateResiduals() {
			for (int i = 0; i < fData.length; ++i)
				fResiduals[i] = fData[i] - fSeasonal[i] - fTrend[i];
		}

		/**
		 * Compute the residual-based weights used in the robustness iterations.
		 */
		private void computeResidualWeights() {

			// TODO: There can be problems if "robust" iterations are done but MAD ~= 0. May want to put a floor on c001.

			// The residual-based weights are a "bisquare" weight based on the residual deviation compared to 6 times the
			// median absolute deviation (MAD). First compute 6 * MAD. (The sort could be a selection but this is
			// not critical as the rest of the algorithm is higher complexity.)

			for (int i = 0; i < fData.length; ++i)
				fWeights[i] = Math.abs(fData[i] - fSeasonal[i] - fTrend[i]);

			Arrays.sort(fWeights);

			// For an even number of elements, the median is the average of the middle two.
			// With proper indexing this formula works either way at the cost of some
			// superfluous work when the number is odd.

			final int mi0 = (fData.length + 1) / 2 - 1; // n = 5, mi0 = 2; n = 4, mi0 = 1
			final int mi1 = fData.length - mi0 - 1;        // n = 5, mi1 = 2; n = 4, mi1 = 2

			final double sixMad = 3.0 * (fWeights[mi0] + fWeights[mi1]);
			final double c999 = 0.999 * sixMad;
			final double c001 = 0.001 * sixMad;

			for (int i = 0; i < fData.length; ++i) {
				double r = Math.abs(fData[i] - fSeasonal[i] - fTrend[i]);
				if (r <= c001) {
					fWeights[i] = 1.0;
				} else if (r <= c999) {
					final double h = r / sixMad;
					final double w = 1.0 - h * h;
					fWeights[i] = w * w;
				} else {
					fWeights[i] = 0.0;
				}
			}
		}

		/**
		 * Smooth the STL seasonal component with quadratic LOESS and recompute the residual.
		 *
		 * @param width             the width of the LOESS smoother used to smooth the seasonal component.
		 */
		public void smoothSeasonal(int width) {
			smoothSeasonal(width, true);
		}

		/**
		 * Smooth the STL seasonal component with quadratic LOESS and recompute the residual.
		 *
		 * @param width             the width of the LOESS smoother used to smooth the seasonal component.
		 * @param restoreEndPoints  whether to restore the endpoints to their original values
		 * @param useWeights        whether to use the STL weights in the seasonal smoothing
		 */
		public void smoothSeasonal(int width, boolean restoreEndPoints) {

			// Ensure that LOESS smoother width is odd and >= 3.

			width = Math.max(3, width);
			if (width % 2 == 0)
				++width;

			// Quadratic smoothing of the seasonal component.
			// Do NOT perform linear interpolation between smoothed points - the quadratic spline can accommodate
			// sharp changes and linear interpolation would cut off peaks/valleys.

			LoessSmoother.Builder builder = new LoessSmoother.Builder().setWidth(width);
			builder.setDegree(2);
			builder.setJump(1);

			LoessSmoother seasonalSmoother = builder.setData(fSeasonal).build();
			double[] smoothedSeasonal = seasonalSmoother.smooth();

			// TODO: Calculate the variance reduction in smoothing the seasonal.

			// Update the seasonal with the smoothed values.

			// TODO: This is not very good - it causes discontinuities a the endpoints.
			//       Better to transition to linear in the last half-smoother width.

			// Restore the end-point values as the smoother will tend to over-modify these.

			double s0 = fSeasonal[0];
			double sN = fSeasonal[fSeasonal.length - 1];
			System.arraycopy(smoothedSeasonal, 0, fSeasonal, 0, smoothedSeasonal.length);

			if (restoreEndPoints) {
				fSeasonal[0] = s0;
				fSeasonal[fSeasonal.length - 1] = sN;
			}

			for (int i = 0; i < smoothedSeasonal.length; ++i)
				fResiduals[i] = fData[i] - fTrend[i] - fSeasonal[i];
		}
	}

	/**
	 * Decompose the input data into the seasonal, trend and residual components.
	 *
	 * @return the STL decomposition
	 */
	public Decomposition decompose() {
		// TODO: Pass input data to decompose and reallocate buffers based on that size.

		fDecomposition = new Decomposition(fData);

		int outerIteration = 0;
		while (true) {

			boolean useResidualWeights = outerIteration > 0;

			for (int iteration = 0; iteration < fInnerIterations; ++iteration) {
				smoothSeasonalSubCycles(useResidualWeights);
				removeSeasonality();
				updateSeasonalAndTrend(useResidualWeights);
			}

			if (++outerIteration > fRobustIterations)
				break;

			fDecomposition.computeResidualWeights();
		}

		fDecomposition.updateResiduals();

		Decomposition result = fDecomposition;

		fDecomposition = null;

		return result;
	}

	/**
	 * The seasonal component is computed by doing smoothing on the cyclic sub-series after removing the trend. The
	 * current estimate of the trend is removed, then the detrended data is separated into sub-series (e.g. all the
	 * Januaries, all the Februaries, etc., for yearly data), and these sub-series are smoothed and extrapolated into
	 * fExtendedSeasonal.
	 */
	private void smoothSeasonalSubCycles(boolean useResidualWeights) {

		double[] data = fDecomposition.fData;
		double[] trend = fDecomposition.fTrend;
		double[] weights = fDecomposition.fWeights;

		for (int i = 0; i < data.length; ++i)
			fDetrend[i] = data[i] - trend[i];

		double[] residualWeights = useResidualWeights ? weights : null;

		fCyclicSubSeriesSmoother.smoothSeasonal(fDetrend, fExtendedSeasonal, residualWeights);
	}

	/**
	 * The lowpass calculation takes the extended seasonal results and smooths them with three moving averages and a
	 * LOESS smoother to remove the seasonality.
	 */
	private void removeSeasonality() {
		// TODO: This creates some garbage - see if its a problem. If so we could preallocate these work arrays and
		// change the code to reuse them.

		// The moving average "erodes" data from the boundaries. We start with:
		//
		// extendedSeasonal.length == data.length + 2 * periodicity
		//
		// and the length after each pass is.................................
		double[] pass1 = simpleMovingAverage(fExtendedSeasonal, fPeriodLength); // data.length + periodLength + 1
		double[] pass2 = simpleMovingAverage(pass1, fPeriodLength);                // data.length + 2
		double[] pass3 = simpleMovingAverage(pass2, 3);                            // data.length

		// assert pass3.length == fData.length; // testing sanity check.

		LoessSmoother lowPassLoess = fLowpassLoessFactory.setData(pass3).build();
		fDeSeasonalized = lowPassLoess.smooth();

		// dumpDebugData("lowpass", fDeSeasonalized);
	}

	/**
	 * The new seasonal component is computed by removing the low-pass smoothed seasonality from the extended
	 * seasonality, and the trend is recalculated by subtracting this new seasonality from the data.
	 */
	private void updateSeasonalAndTrend(boolean useResidualWeights) {

		double[] data = fDecomposition.fData;
		double[] trend = fDecomposition.fTrend;
		double[] weights = fDecomposition.fWeights;
		double[] seasonal = fDecomposition.fSeasonal;

		for (int i = 0; i < data.length; ++i) {
			seasonal[i] = fExtendedSeasonal[fPeriodLength + i] - fDeSeasonalized[i];
			trend[i] = data[i] - seasonal[i];
		}

		// dumpDebugData("seasonal", seasonal);
		// dumpDebugData("trend0", trend);

		double[] residualWeights = useResidualWeights ? weights : null;

		LoessSmoother trendSmoother = fLoessSmootherFactory.setData(trend).setExternalWeights(residualWeights).build();

		System.arraycopy(trendSmoother.smooth(), 0, trend, 0, trend.length);

		// dumpDebugData("trend", trend);
	}

	@Override
	public String toString() {
		return String.format(
				"SeasonalTrendLoess: [\n" +
						"inner iterations     = %d\n" +
						"outer iterations     = %d\n" +
						"periodicity          = %d\n" +
						"seasonality settings = %s\n" +
						"trend settings       = %s\n" +
						"lowpass settings     = %s\n]",
				this.fInnerIterations, this.fRobustIterations, this.fPeriodLength,
				this.fSeasonalSettings, this.fTrendSettings, this.fLowpassSettings);
	}

//	@SuppressWarnings("unused")
//	private static void dumpDebugData(String name, double[] data) {
//		for (int i = 0; i < data.length; ++i) {
//			System.out.println(String.format("%s(%3d) = %22.15e", name, i, data[i]));
//		}
//	}
}

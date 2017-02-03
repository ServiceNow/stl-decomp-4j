package com.snc.ds.stats.stl;

import static com.snc.ds.stats.TimeSeriesUtilities.simpleMovingAverage;

import java.util.Arrays;

/**
 * Java implementation of the Seasonal-Trend-Loess algorithm for evenly spaced data. This is basically a direct port of
 * the RATFOR from the netlib stl package.
 *
 * Created by Jim Crotinger on 18-Apr-2016.
 */
public class SeasonalTrendLoess {

	private final double[] fData;
	private final int fPeriodLength;
	private final LoessSettings fSeasonalSettings;
	private final LoessSettings fTrendSettings;
	private final LoessSettings fLowpassSettings;
	private final int fInnerIterations;
	private final int fRobustIterations;

	private final double[] fTrend;
	private final double[] fSeasonal;
	private final double[] fResiduals;
	private final double[] fWeights;

	private final double[] fDetrend;
	private final double[] fExtendedSeasonal;

	private double[] fDeSeasonalized; // TODO: Garbage - can this be made in-place?

	private final CyclicSubSeriesSmoother cyclicSubSeries;

	private static int calcDefaultTrendWidth(int periodicity, int seasonalWidth) {
		// This formula is based on a numerical stability analysis in the original paper.
		return (int) (1.5 * periodicity / (1 - 1.5 / seasonalWidth) + 0.5);
	}

	/**
	 * Construct STL specifying full details of the LOESS smoothers via LoessSettings objects.
	 *
	 * @param data
	 *            - the data to be decomposed
	 * @param periodicity
	 *            - the periodicity of the data
	 * @param ni
	 *            - the number of inner iterations
	 * @param no
	 *            - the number of outer "robustness" iterations
	 * @param seasonalSettings
	 *            - the settings for the LOESS smoother for the cyclic sub-series
	 * @param trendSettings
	 *            - the settings for the LOESS smoother for the trend component
	 * @param lowpassSettings
	 *            - the settings for the LOESS smoother used in de-seasonalizing
	 */
	public SeasonalTrendLoess(double[] data, int periodicity, int ni, int no, LoessSettings seasonalSettings,
			LoessSettings trendSettings, LoessSettings lowpassSettings) {

		if (periodicity < 2)
			throw new RuntimeException("periodicity must be at least 2");

		if (data.length < 2 * periodicity)
			throw new RuntimeException("Data series must be at least 2 * periodicity in length");

		this.fData = data;
		this.fPeriodLength = periodicity;
		this.fSeasonalSettings = seasonalSettings;
		this.fTrendSettings = trendSettings;
		this.fLowpassSettings = lowpassSettings;
		this.fInnerIterations = ni;
		this.fRobustIterations = no;
		this.cyclicSubSeries = new CyclicSubSeriesSmoother(seasonalSettings, data.length, periodicity);

		final int size = fData.length;
		fTrend = new double[size];
		fSeasonal = new double[size];
		fResiduals = new double[size];
		fWeights = new double[size];
		fDetrend = new double[size];
		fExtendedSeasonal = new double[size + 2 * fPeriodLength];

		Arrays.fill(fWeights, 1.0);
	}

	/**
	 * Construct STL specifying widths and iterations explicitly.
	 *
	 * @param data
	 * @param periodicity
	 * @param seasonalWidth
	 * @param trendWidth
	 * @param ni
	 * @param no
	 */
	@SuppressWarnings("WeakerAccess")
	public SeasonalTrendLoess(double[] data, int periodicity, int seasonalWidth, int trendWidth, int ni, int no) {
		this(data, periodicity, ni, no, new LoessSettings(seasonalWidth), new LoessSettings(trendWidth),
				new LoessSettings(periodicity));
	}

	/**
	 * Construct STL specifying seasonal width and iterations. Trend width calculated.
	 *
	 * @param data
	 * @param periodicity
	 * @param seasonalWidth
	 * @param ni
	 * @param no
	 */
	public SeasonalTrendLoess(double[] data, int periodicity, int seasonalWidth, int ni, int no) {
		this(data, periodicity, seasonalWidth, calcDefaultTrendWidth(periodicity, seasonalWidth), ni, no);
	}

	/**
	 * Construct STL specifying full settings for the seasonal LOESS, and iterations. Trend width calculated.
	 *
	 * @param data
	 * @param periodicity
	 * @param seasonalSettings
	 * @param ni
	 * @param no
	 */
	@SuppressWarnings("WeakerAccess")
	public SeasonalTrendLoess(double[] data, int periodicity, LoessSettings seasonalSettings, int ni, int no) {
		this(data, periodicity,
				ni, no,
				seasonalSettings,
				new LoessSettings(calcDefaultTrendWidth(periodicity, seasonalSettings.getWidth())),
				new LoessSettings(periodicity));
	}

	/**
	 * Construct STL specifying periodicity and smoothing widths, with iterations chosen via the robust flag.
	 *
	 * @param data
	 * @param periodicity
	 * @param seasonalWidth
	 * @param trendWidth
	 * @param robust
	 */
	@SuppressWarnings("WeakerAccess")
	public SeasonalTrendLoess(double[] data, int periodicity, int seasonalWidth, int trendWidth, boolean robust) {
		this(data, periodicity, seasonalWidth, trendWidth, robust ? 1 : 2, robust ? 15 : 0);
	}

	/**
	 * Construct STL specifying periodicity and seasonal width. Computes trend-smoothing width. Iterations chosen via
	 * robust flag.
	 *
	 * @param data
	 * @param periodicity
	 * @param seasonalWidth
	 * @param robust
	 */
	public SeasonalTrendLoess(double[] data, int periodicity, int seasonalWidth, boolean robust) {
		this(data, periodicity, seasonalWidth, calcDefaultTrendWidth(periodicity, seasonalWidth), robust);
	}

	/**
	 * Factory method to perform a non-robust STL decomposition enforcing strict periodicity.
	 * <p/>
	 * Meant for diagnostic purposes only.
	 *
	 * @param data
	 *            the data to analyze
	 * @param periodicity
	 *            the (suspected) periodicity of the data
	 * @return SeasonalTrendLoess object with the decomposition already performed.
	 */
	public static SeasonalTrendLoess performPeriodicDecomposition(
			double[] data,
			@SuppressWarnings("SameParameterValue") int periodicity
	) {
		// The LOESS interpolator with degree 0 and a very long window (arbitrarily chosen to be 100 times the length of
		// the array) will interpolate all points as the average value of the series. This particular setting is used
		// for smoothing the seasonal sub-cycles, so the end result is that the seasonal component of the decomposition
		// is exactly periodic.
		int width = 100 * data.length;
		int degree = 0;
		LoessSettings seasonalSettings = new LoessSettings(width, degree);

		// This fit is for diagnostic purposes, so we just do a single inner iteration.
		int ni = 1;
		int no = 0;
		SeasonalTrendLoess stl = new SeasonalTrendLoess(data, periodicity, seasonalSettings, ni, no);
		stl.decompose();
		return stl;
	}

	/**
	 * Factory method to perform a (somewhat) robust STL decomposition enforcing strict periodicity.
	 * <p/>
	 * Meant for diagnostic purposes only.
	 *
	 * @param data
	 *            the data to analyze
	 * @param periodicity
	 *            the (suspected) periodicity of the data
	 * @return SeasonalTrendLoess object with the decomposition already performed.
	 */
	public static SeasonalTrendLoess performRobustPeriodicDecomposition(
			double[] data,
			@SuppressWarnings("SameParameterValue") int periodicity
	) {

		// The LOESS interpolator with degree 0 and a very long window (arbitrarily chosen to be 100 times the length of
		// the array) will interpolate all points as the average value of the series. This particular setting is used
		// for smoothing the seasonal sub-cycles, so the end result is that the seasonal component of the decomposition
		// is exactly periodic.
		int width = 100 * data.length;
		int degree = 0;
		LoessSettings seasonalSettings = new LoessSettings(width, degree);

		// This fit is for diagnostic purposes, so we just do a single outer iteration.
		int ni = 1;
		int no = 1;
		SeasonalTrendLoess stl = new SeasonalTrendLoess(data, periodicity, seasonalSettings, ni, no);
		stl.decompose();
		return stl;
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
	public double[] getResiduals() {
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

	/**
	 * Decompose the input data into the seasonal, trend and residual components.
	 */
	public void decompose() {
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

			computeResidualWeights();
		}

		for (int i = 0; i < fData.length; ++i)
			fResiduals[i] = fData[i] - fSeasonal[i] - fTrend[i];
	}

	/**
	 * Smooth the STL seasonal component with quadratic LOESS and recompute the residual.
	 *
	 * @param width
	 *            the width of the LOESS smoother used to smooth the seasonal component.
	 */
	public void smoothSeasonal(int width) {
		int degree = 2; // quadratic
		// Don't use linear interpolation - the quadratic spline can accommodate sharp changes and the linear
		// interpolation will cut off peaks/valleys.
		int jump = 1;
		LoessSettings settings = new LoessSettings(width, degree, jump);
		LoessSmoother seasonalSmoother = new LoessSmoother(settings, fSeasonal, null);
		double[] smoothedSeasonal = seasonalSmoother.smooth();

		// Update the seasonal with the smoothed values.

		// Restore the end-point values as the smoother will tend to over-modify these.

		double s0 = fSeasonal[0];
		double sN = fSeasonal[fSeasonal.length - 1];
		System.arraycopy(smoothedSeasonal, 0, fSeasonal, 0, smoothedSeasonal.length);

		fSeasonal[0] = s0;
		fSeasonal[fSeasonal.length - 1] = sN;

		for (int i = 0; i < smoothedSeasonal.length; ++i)
			fResiduals[i] = fData[i] - fTrend[i] - fSeasonal[i];

	}

	/**
	 * The seasonal component is computed by doing smoothing on the cyclic sub-series after removing the trend. The
	 * current estimate of the trend is removed, then the detrended data is separated into sub-series (e.g. all the
	 * Januaries, all the Februaries, etc., for yearly data), and these sub-series are smoothed and extrapolated into
	 * fExtendedSeasonal.
	 */
	private void smoothSeasonalSubCycles(boolean useResidualWeights) {

		for (int i = 0; i < fData.length; ++i)
			fDetrend[i] = fData[i] - fTrend[i];

		double[] residualWeights = useResidualWeights ? fWeights : null;

		cyclicSubSeries.smoothSeasonal(fDetrend, fExtendedSeasonal, residualWeights);
	}

	/**
	 * The lowpass calculation takes the extended seasonal results and smoothes them with three moving averages and a
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
		double[] pass2 = simpleMovingAverage(pass1, fPeriodLength);				// data.length + 2
		double[] pass3 = simpleMovingAverage(pass2, 3);							// data.length

		// assert pass3.length == fData.length; // testing sanity check.

		LoessSmoother lowPassLoess = new LoessSmoother(fLowpassSettings, pass3, null);
		fDeSeasonalized = lowPassLoess.smooth();

		// dumpDebugData("lowpass", lowpass);
	}

	/**
	 * The new seasonal component is computed by removing the low-pass smoothed seasonality from the extended
	 * seasonality, and the trend is recalculated by subtracting this new seasonality from the data.
	 */
	private void updateSeasonalAndTrend(boolean useResidualWeights) {

		for (int i = 0; i < fData.length; ++i) {
			fSeasonal[i] = fExtendedSeasonal[fPeriodLength + i] - fDeSeasonalized[i];
			fTrend[i] = fData[i] - fSeasonal[i];
		}

		// dumpDebugData("seasonal", seasonal);
		// dumpDebugData("trend0", trend);

		double[] residualWeights = useResidualWeights ? fWeights : null;
		LoessSmoother trendSmoother = new LoessSmoother(fTrendSettings, fTrend, residualWeights);

		System.arraycopy(trendSmoother.smooth(), 0, fTrend, 0, fTrend.length);

		// dumpDebugData("trend", trend);
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
		final int mi1 = fData.length - mi0 - 1;		// n = 5, mi1 = 2; n = 4, mi1 = 2

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

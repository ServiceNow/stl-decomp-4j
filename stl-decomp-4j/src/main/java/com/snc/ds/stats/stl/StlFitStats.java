package com.snc.ds.stats.stl;

/**
 * StlFitStats analyzes an STL decomposition, computing basic statistics for the original data and the decomposition.
 *
 * Created by Jim Crotinger on 1-Jun-2016
 */
@SuppressWarnings("WeakerAccess")
public class StlFitStats {
	private final int fSampleSize;
	private final double fDataMean;
	private final double fDataVariance;
	private final double fTrendMean;
	private final double fTrendRange;
	private final double fSeasonalMean;
	private final double fSeasonalVariance;
	private final double fResidualMean;
	private final double fResidualVariance;
	private final double fDeSeasonalMean;
	private final double fDeSeasonalVariance;
	private final double fDeTrendMean;
	private final double fDeTrendVariance;
	private final double fSeasonalRange;

	/**
	 * Analyze the STL decomposition, computing basic statistics for the original data and the decomposition.
	 *
	 * @param stl a SeasonalTrendLoess object after the decomposition is performed
	 */
	public StlFitStats(SeasonalTrendLoess stl) {
		int length = stl.getData().length;
		// Unnecessary since STL guarantees this, so it can't be tested:
		// Preconditions.checkArgument(length >= 4, "STL Decomposition must have at least 4 data points");

		final double[] data = stl.getData();
		final double[] trend = stl.getTrend();
		final double[] seasonal = stl.getSeasonal();
		final double[] residuals = stl.getResiduals();

		double dataSum = 0.0;
		double dataSqSum = 0.0;
		double trendSum = 0.0;
		double trendMax = -1.0e100;
		double trendMin = 1.0e100;
		double seasonalSum = 0.0;
		double seasonalSqSum = 0.0;
		double seasonalMax = -1.0e100;
		double seasonalMin = 1.0e100;
		double residualSum = 0.0;
		double residualSqSum = 0.0;
		double deSeasonalSum = 0.0;
		double deSeasonalSqSum = 0.0;
		double deTrendSum = 0.0;
		double deTrendSqSum = 0.0;

		for (int i = 0; i < length; ++i) {
			double d = data[i];
			double t = trend[i];
			double s = seasonal[i];
			double r = residuals[i];
			double f = d - s;
			double dt = d - t;

			dataSum += d;
			dataSqSum += d * d;

			trendSum += t;
			if (t > trendMax)
				trendMax = t;
			if (t < trendMin)
				trendMin = t;

			seasonalSum += s;
			seasonalSqSum += s * s;
			if (s > seasonalMax)
				seasonalMax = s;
			if (s < seasonalMin)
				seasonalMin = s;

			residualSum += r;
			residualSqSum += r * r;

			deSeasonalSum += f;
			deSeasonalSqSum += f * f;

			deTrendSum += dt;
			deTrendSqSum += dt * dt;
		}

		double denom = 1.0 / length;

		fDataMean = dataSum * denom;
		fTrendMean = trendSum * denom;
		fSeasonalMean = seasonalSum * denom;
		fResidualMean = residualSum * denom;
		fDeSeasonalMean = deSeasonalSum * denom;
		fDeTrendMean = deTrendSum * denom;

		// The data is from a valid STL decomposition, so length = 4 at minimum.

		double corrBC = length / (length - 1.0); // Bessel's correction
		double denomBC = 1.0 / (length - 1.0);
		fDataVariance = dataSqSum * denomBC - fDataMean * fDataMean * corrBC;
		fTrendRange = trendMax - trendMin;
		fSeasonalVariance = seasonalSqSum * denomBC - fSeasonalMean * fSeasonalMean * corrBC;
		fSeasonalRange = seasonalMax - seasonalMin;
		fResidualVariance = residualSqSum * denomBC - fResidualMean * fResidualMean * corrBC;
		fDeSeasonalVariance = deSeasonalSqSum * denomBC - fDeSeasonalMean * fDeSeasonalMean * corrBC;
		fDeTrendVariance = deTrendSqSum * denomBC - fDeTrendMean * fDeTrendMean * corrBC;

		fSampleSize = length;
	}

	/**
	 * Get the mean of the trend.
	 *
	 * @return the mean value of the trend
	 */
	public double getTrendMean() {
		return fTrendMean;
	}

	/**
	 * Get the range (max - min) of the trend.
	 *
	 * @return the range of the trend
	 */
	public double getTrendRange() {
		return fTrendRange;
	}

	/**
	 * Get the mean of the data.
	 *
	 * @return the mean value of the data
	 */
	public double getDataMean() {
		return fDataMean;
	}

	/**
	 * Get the variance of the data.
	 *
	 * @return the variance of the data
	 */
	public double getDataVariance() {
		return fDataVariance;
	}

	/**
	 * Get the standard deviation of the data.
	 *
	 * @return the standard deviation of the data
	 */
	public double getDataStdDev() {
		return Math.sqrt(fDataVariance);
	}

	/**
	 * Get the mean of the seasonal component. Should be near zero if the length is an even multiple of the period.
	 *
	 * @return the mean value of the seasonal component
	 */
	public double getSeasonalMean() {
		return fSeasonalMean;
	}

	/**
	 * Get the variance of the seasonal component.
	 *
	 * @return the variance of the seasonal component
	 */
	public double getSeasonalVariance() {
		return fSeasonalVariance;
	}

	/**
	 * Get the standard deviation of the seasonal component.
	 *
	 * @return the standard deviation of the seasonal component
	 */
	public double getSeasonalStdDev() {
		return Math.sqrt(fSeasonalVariance);
	}

	/**
	 * Get the range (max - min) of the seasonal component.
	 *
	 * @return the range of the seasonal component
	 */
	public double getSeasonalRange() {
		return fSeasonalRange;
	}

	/**
	 * Get the mean of the residual. Should be near zero.
	 *
	 * @return the mean of the residual
	 */
	public double getResidualMean() {
		return fResidualMean;
	}

	/**
	 * Get the variance of the residual.
	 *
	 * @return the variance of the residual
	 */
	public double getResidualVariance() {
		return fResidualVariance;
	}

	/**
	 * Get the standard deviation of the residual.
	 *
	 * @return the standard deviation of the residual
	 */
	public double getResidualStdDev() {
		return Math.sqrt(fResidualVariance);
	}

	/**
	 * Get the deseasonalized mean.
	 *
	 * @return the mean of the data after the seasonal component is removed.
	 */
	public double getDeSeasonalMean() {
		return fDeSeasonalMean;
	}

	/**
	 * Get the deseasonalized variance. This is the same as the trend-free variance since if the trend is constant it
	 * will drop out of the variance.
	 *
	 * @return the variance of the data after the seasonal component is removed.
	 */
	public double getDeSeasonalVariance() {
		return fDeSeasonalVariance;
	}


	/**
	 * Get the deseasonalized mean.
	 *
	 * @return the mean of the data after the seasonal component is removed.
	 */
	public double getDeTrendMean() {
		return fDeTrendMean;
	}

	/**
	 * Get the deseasonalized variance. This is the same as the trend-free variance since if the trend is constant it
	 * will drop out of the variance.
	 *
	 * @return the variance of the data after the seasonal component is removed.
	 */
	public double getDeTrendVariance() {
		return fDeTrendVariance;
	}

	/**
	 * Get the estimated variance of the residual sample variance. Based on the assumption that the residuals are
	 * normal.
	 *
	 * @return the estimate of the variance of the residual sample variance.
	 */
	public double getEstimatedVarianceOfResidualSampleVariance() {
		double v = getResidualVariance();
		return 2 * v * v / (fSampleSize - 1);
	}

	/**
	 * Get a Z-Score for the residual variance with no trend relative to the statistics from the STL residual.
	 */
	public double getTrendinessZScore() {
		double resVarVar = getEstimatedVarianceOfResidualSampleVariance();

		return (fDeSeasonalVariance - fResidualVariance) / Math.sqrt(Math.max(1.0e-12, resVarVar));
	}

	/**
	 * Get a Z-Score for the residual variance with no seasonality relative to the statistics from the STL residual.
	 */
	public double getSeasonalZScore() {
		double resVarVar = getEstimatedVarianceOfResidualSampleVariance();

		return (fDeTrendVariance - fResidualVariance) / Math.sqrt(Math.max(1.0e-12, resVarVar));
	}

	@Override
	public String toString() {
		return String.format(
				"Data Mean            = %10f\n" +
				"Data Variance        = %10f\n" +
				"Trend Mean           = %10f\n" +
				"Trend Range          = %10f\n" +
				"Seasonal Mean        = %10f\n" +
				"Seasonal Variance    = %10f\n" +
				"Seasonal Range       = %10f\n" +
				"De-Seasonal Mean     = %10f\n" +
				"De-Seasonal Variance = %10f\n" +
				"De-Trend Mean        = %10f\n" +
				"De-Trend Variance    = %10f\n" +
				"Residual Mean        = %10f\n" +
				"Residual Variance    = %10f\n" +
				"Var(ResSampleVar)    = %10f\n" +
				"Trend Test ZScore    = %10f\n" +
				"Seasonal Test ZScore = %10f\n" +
				"SeasonalVar/ResidVar = %10f",
				getDataMean(), getDataVariance(),
				getTrendMean(), getTrendRange(),
				getSeasonalMean(), getSeasonalVariance(), getSeasonalRange(),
				getDeSeasonalMean(), getDeSeasonalVariance(),
				getDeTrendMean(), getDeTrendVariance(),
				getResidualMean(), getResidualVariance(),
				getEstimatedVarianceOfResidualSampleVariance(),
				getTrendinessZScore(), getSeasonalZScore(),
				getSeasonalVariance() / getResidualVariance());
	}

}

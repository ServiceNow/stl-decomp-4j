package com.snc.ds.stats.stl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Test;

import com.snc.ds.stats.stl.SeasonalTrendLoess.Builder;
import com.snc.ds.stats.stl.SeasonalTrendLoess.Decomposition;

/**
 * Unit test for StlFitStats
 *
 * Created by Jim Crotinger on 1-Jun-2016
 */
public class StlFitStatsTest {

	private final StlTestDataGenerator testDataGenerator = new StlTestDataGenerator();

	@Test
	public void StlStatsSanityTest() {

		int periodicity = 168;
		Builder builder = new Builder().setPeriodLength(periodicity).setSeasonalWidth(2001);
		builder.setInnerIterations(1).setRobustnessIterations(15);

		SeasonalTrendLoess smoother = builder.buildSmoother(testDataGenerator.values);

		Decomposition stl = smoother.decompose();

		StlFitStats stats = new StlFitStats(stl);

		int length = testDataGenerator.values.length;
		SummaryStatistics dataStats = new SummaryStatistics();
		SummaryStatistics trendStats = new SummaryStatistics();
		SummaryStatistics seasonalStats = new SummaryStatistics();
		SummaryStatistics residualStats = new SummaryStatistics();
		SummaryStatistics deSeasonalStats = new SummaryStatistics();

		for (int i = 0; i < length; ++i) {
			dataStats.addValue(stl.getData()[i]);
			trendStats.addValue(stl.getTrend()[i]);
			seasonalStats.addValue(stl.getSeasonal()[i]);
			residualStats.addValue(stl.getResiduals()[i]);
			deSeasonalStats.addValue(stl.getData()[i] - stl.getSeasonal()[i]);
		}

		String statsData = stats.toString();

		assertEquals("raw data mean", dataStats.getMean(), stats.getDataMean(), 1.0e-11);
		assertEquals("raw data variance", dataStats.getVariance(), stats.getDataVariance(), 1.0e-11);
		assertEquals("raw data std dev", dataStats.getStandardDeviation(), stats.getDataStdDev(), 1.0e-11);

		assertEquals("trend mean", trendStats.getMean(), stats.getTrendMean(), 1.0e-11);
		assertEquals("trend range", trendStats.getMax() - trendStats.getMin(), stats.getTrendRange(), 1.0e-11);

		assertEquals("seasonal mean", seasonalStats.getMean(), stats.getSeasonalMean(), 1.0e-11);
		assertEquals("seasonal variance", seasonalStats.getVariance(), stats.getSeasonalVariance(), 1.0e-11);
		assertEquals("seasonal std dev", seasonalStats.getStandardDeviation(), stats.getSeasonalStdDev(), 1.0e-11);

		assertEquals("residual mean", residualStats.getMean(), stats.getResidualMean(), 1.0e-11);
		assertEquals("residual variance", residualStats.getVariance(), stats.getResidualVariance(), 1.0e-11);
		assertEquals("residual std dev", residualStats.getStandardDeviation(), stats.getResidualStdDev(), 1.0e-11);

		assertEquals("deseasonal mean", deSeasonalStats.getMean(), stats.getDeSeasonalMean(), 1.0e-11);
		assertEquals("deseasonal variance", deSeasonalStats.getVariance(), stats.getDeSeasonalVariance(), 1.0e-11);

		double residualVar = stats.getResidualVariance();
		double resSampleVarVar = residualVar * residualVar * 2 / (length - 1);
		double resSampleVarStd = Math.sqrt(resSampleVarVar);
		double trendZ = (stats.getDeSeasonalVariance() - residualVar) / resSampleVarStd;
		assertEquals("trend test z-score", trendZ, stats.getTrendinessZScore(), 1.0e-11);

		assertEquals("Data Mean            =   7.963550\n" + "Data Variance        =  53.996027\n"
				+ "Trend Mean           =   7.691074\n" + "Trend Range          =   1.436350\n"
				+ "Seasonal Mean        =   0.106905\n" + "Seasonal Variance    =  37.427790\n"
				+ "Seasonal Range       =  26.871361\n" + "De-Seasonal Mean     =   7.856645\n"
				+ "De-Seasonal Variance =  16.398024\n" + "De-Trend Mean        =   0.272476\n"
				+ "De-Trend Variance    =  53.741629\n" + "Residual Mean        =   0.165571\n"
				+ "Residual Variance    =  16.161776\n" + "Var(ResSampleVar)    =   0.357079\n"
				+ "Trend Test ZScore    =   0.395354\n" + "Seasonal Test ZScore =  62.888777\n"
				+ "SeasonalVar/ResidVar =   2.315822", statsData);
	}

	@Test
	public void pureSineTest() {

		long seed = 1234567L; // System.nanoTime() // change this to do random stress test

		double[] data = testDataGenerator.createNoisySeasonalData(144, 12, 1.0, 0.0, 0.0, seed); // sin(2 \pi i / 12)

		final Builder builder = new Builder().setPeriodLength(12).setSeasonalWidth(7).setNonRobust();
		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		Decomposition stl = smoother.decompose();

		StlFitStats stats = new StlFitStats(stl);

		SummaryStatistics dataStats = new SummaryStatistics();
		for (double v : stl.getData())
			dataStats.addValue(v);

		assertEquals("raw data mean", 0.0, stats.getDataMean(), 1.0e-11);
		assertEquals("raw data variance", dataStats.getVariance(), stats.getDataVariance(), 1.0e-11);
		assertEquals("raw data std dev", dataStats.getStandardDeviation(), stats.getDataStdDev(), 1.0e-11);

		assertEquals("trend mean", 0.0, stats.getTrendMean(), 1.0e-11);
		assertEquals("trend range", 0.0, stats.getTrendRange(), 1.0e-11);

		assertEquals("seasonal mean", 0.0, stats.getSeasonalMean(), 1.0e-11);
		assertEquals("seasonal variance", dataStats.getVariance(), stats.getSeasonalVariance(), 1.0e-11);
		assertEquals("seasonal std dev", dataStats.getStandardDeviation(), stats.getSeasonalStdDev(), 1.0e-11);

		assertEquals("residual mean", 0.0, stats.getResidualMean(), 1.0e-11);
		assertEquals("residual variance", 0.0, stats.getResidualVariance(), 1.0e-11);
		assertEquals("residual std dev", 0.0, stats.getResidualStdDev(), 1.0e-11);

		assertEquals("deseasonal mean", 0.0, stats.getDeSeasonalMean(), 1.0e-11);
		assertEquals("deseasonal variance", 0.0, stats.getDeSeasonalVariance(), 1.0e-11);

		assertEquals("trend test z-score", 0.0, stats.getTrendinessZScore(), 1.0e-11);
		assertEquals("seasonal test z-score", stats.getSeasonalVariance(), 1.0e-6 * stats.getSeasonalZScore(), 1.0e-11);
	}

	@Test
	public void pureTrendTest() {

		long seed = 1234567L; // System.nanoTime() // change this to do random stress test

		// linear trend (2 \pi i / 12)
		double[] data = testDataGenerator.createNoisySeasonalData(144, 12, 0.0, 1.0, 0.0, seed);

		final Builder builder = new Builder().setPeriodLength(12).setSeasonalWidth(7).setNonRobust();
		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		Decomposition stl = smoother.decompose();

		StlFitStats stats = new StlFitStats(stl);

		SummaryStatistics dataStats = new SummaryStatistics();
		for (double v : stl.getData())
			dataStats.addValue(v);

		assertEquals("raw data mean", dataStats.getMean(), stats.getDataMean(), 1.0e-11);
		assertEquals("raw data variance", dataStats.getVariance(), stats.getDataVariance(), 1.0e-11);
		assertEquals("raw data std dev", dataStats.getStandardDeviation(), stats.getDataStdDev(), 1.0e-11);

		assertEquals("trend mean", dataStats.getMean(), stats.getTrendMean(), 1.0e-11);
		assertEquals("trend range", 2.0 * dataStats.getMean(), stats.getTrendRange(), 1.0e-11);

		assertEquals("seasonal mean", 0.0, stats.getSeasonalMean(), 1.0e-11);
		assertEquals("seasonal variance", 0.0, stats.getSeasonalVariance(), 1.0e-11);
		assertEquals("seasonal std dev", 0.0, stats.getSeasonalStdDev(), 1.0e-11);

		assertEquals("residual mean", 0.0, stats.getResidualMean(), 1.0e-11);
		assertEquals("residual variance", 0.0, stats.getResidualVariance(), 1.0e-11);
		assertEquals("residual std dev", 0.0, stats.getResidualStdDev(), 1.0e-11);

		assertEquals("deseasonal mean", dataStats.getMean(), stats.getDeSeasonalMean(), 1.0e-11);
		assertEquals("deseasonal variance", dataStats.getVariance(), stats.getDeSeasonalVariance(), 1.0e-11);
		assertEquals("trend test z-score", dataStats.getVariance(), 1.0e-6 * stats.getTrendinessZScore(), 1.0e-11);
		assertEquals("seasonal test z-score", 0.0, 1.0e-6 * stats.getSeasonalZScore(), 1.0e-11);
	}

	@Test
	public void noisyTrendyTest() {
		double maxSmoothedZScore = 0;
		double maxFracVar = 0;
		int trials = 1000;

		long seed = 1234567L; // System.nanoTime() // change this to do random stress test

		for (int i = 0; i < trials; ++i) {
			double[] data = testDataGenerator.createNoisySeasonalData(168 * 4, 168, 0.0, 0.2, 1.0, seed++);

			Decomposition stl = SeasonalTrendLoess.performPeriodicDecomposition(data, 168);

			StlFitStats stats = new StlFitStats(stl);

			assertTrue(String.format("iteration %d", i), stats.getTrendinessZScore() > 3.0);

			stl.smoothSeasonal(15);
			StlFitStats smoothedStats = new StlFitStats(stl);

			double fractionalVariance = smoothedStats.getSeasonalVariance() / smoothedStats.getDeTrendVariance();
			if (fractionalVariance > maxFracVar)
				maxFracVar = fractionalVariance;

			assertTrue(
					String.format("Seasonal Variance Check: seed %d; iteration %d; stats = %s", seed, i, smoothedStats),
					fractionalVariance < 0.33);

			double seasonalZScore = smoothedStats.getSeasonalZScore();
			if (seasonalZScore > maxSmoothedZScore)
				maxSmoothedZScore = seasonalZScore;

			assertTrue(
					String.format("Seasonal Z-Score Check: seed %d; iteration %d; stats = %s", seed, i, smoothedStats),
					seasonalZScore < 3.0);
		}
	}

	@Test
	public void noisySeasonalTest() {
		long seed = 1234567L; // System.nanoTime() // change this to do random stress test
		int numAverages = 1; // Increase to tighten the statistics on the sample statistics
		int trials = 100;
		double start = 1.5;
		double delta = 0.0;

		SummaryStatistics seasonalZScoreStats = new SummaryStatistics();
		SummaryStatistics varianceFractionStats = new SummaryStatistics();

		SummaryStatistics fractionClassifiedStats = new SummaryStatistics();
		SummaryStatistics sampleZScoreStats = new SummaryStatistics();
		SummaryStatistics sampleVarFracStats = new SummaryStatistics();

		for (int j = 0; j < numAverages; ++j) {
			int count = 0;
			double seasonalAmplitude = start + delta * j;
			double noiseSigma = 3.0;

			for (int i = 0; i < trials; ++i) {

				double[] data = testDataGenerator.createNoisySeasonalData(
						168 * 4, 168, seasonalAmplitude, 0.0, noiseSigma, seed++);

				Decomposition stl = SeasonalTrendLoess.performRobustPeriodicDecomposition(data, 168);

				StlFitStats stats = new StlFitStats(stl);

				assertTrue(stats.getTrendinessZScore() < 3.0);

				stl.smoothSeasonal(15);
				StlFitStats smoothedStats = new StlFitStats(stl);

				double vf = smoothedStats.getSeasonalVariance() / smoothedStats.getResidualVariance();
				varianceFractionStats.addValue(vf);

				double z = smoothedStats.getSeasonalZScore();
				seasonalZScoreStats.addValue(z);

				if (z > 3.0)
					++count;
			}

			double rate = ((double) count) / trials;

			fractionClassifiedStats.addValue(rate);

			sampleZScoreStats.addValue(seasonalZScoreStats.getMean());
			sampleVarFracStats.addValue(varianceFractionStats.getMean());

			seasonalZScoreStats.clear();
			varianceFractionStats.clear();
		}

		assertTrue("Min Average Z-Score", sampleZScoreStats.getMin() > 3.13);
		assertTrue("Avg Average Z-Score", Math.abs(sampleZScoreStats.getMean() - 3.64) < 0.06);
		assertTrue("Max Average Z-Score", sampleZScoreStats.getMax() < 4.13);

		assertTrue("Min Var Frac", sampleVarFracStats.getMin() > 0.173);
		assertTrue("Avg Var Frac", Math.abs(sampleVarFracStats.getMean() - 0.193) < 0.01);
		assertTrue("Max Var Frac", sampleVarFracStats.getMax() < 0.213);
	}

	@Test
	public void noisyTrendyBadSeedTest() throws IOException {

		long seed = 16951029831410L; // 18823810773670L

		double[] data = testDataGenerator.createNoisySeasonalData(168 * 4, 168, 0.0, 0.2, 1.0, 16951029831410L);

		Decomposition stl = SeasonalTrendLoess.performRobustPeriodicDecomposition(data, 168);

		StlFitStats stats = new StlFitStats(stl);

		assertTrue(stats.getTrendinessZScore() > 3.0);

		stl.smoothSeasonal(15);
		StlFitStats smoothedStats = new StlFitStats(stl);

		assertTrue(String.format("Seasonal Variance Check: seed %d; stats = %s", seed, smoothedStats),
				smoothedStats.getSeasonalVariance() < 0.12 * smoothedStats.getDeTrendVariance());
		assertTrue(String.format("Seasonal Z-Score Check: seed %d; stats = %s", seed, smoothedStats),
				smoothedStats.getSeasonalZScore() < 3.0); // TODO: A better seasonal Z-Score
	}

}

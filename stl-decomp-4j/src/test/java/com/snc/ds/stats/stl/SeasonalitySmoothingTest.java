package com.snc.ds.stats.stl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Ignore;
import org.junit.Test;

public class SeasonalitySmoothingTest {

	private final StlTestDataGenerator fTestData = new StlTestDataGenerator();

	@Test
	public void smoothSeasonalityTest() throws IOException {

		final double[] data = fTestData.values;

		int periodicity = 168;
		SeasonalTrendLoess stl = new SeasonalTrendLoess(data, periodicity, 2001, 1, 15);

		stl.decompose();

		double epsilon = 5.0e-03; // Doing absolute differences and some values > 10

		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResiduals();

		for (int i = 0; i < data.length; ++i) {
			assertEquals(String.format("seasonal[%d]", i), fTestData.trend[i], trend[i], epsilon);
			assertEquals(String.format("trend[%d]", i), fTestData.seasonal[i], seasonal[i], 20 * epsilon);
			assertEquals(String.format("residuals[%d]", i), fTestData.residuals[i], residuals[i],
					20 * epsilon);
		}

		double[] extendedSeasonal = new double[seasonal.length + 2 * periodicity];
		LoessSettings seasonalSettings = new LoessSettings(2001);
		CyclicSubSeriesSmoother seasonal_smoother = new CyclicSubSeriesSmoother(seasonalSettings, seasonal.length, periodicity);
		seasonal_smoother.smoothSeasonal(seasonal, extendedSeasonal, null);

		double[] modelSeasonal = new double[2*periodicity];
//		double[] modelTrend = new double[2*periodicity];

//		double trendEnd = trend[trend.length - 1];
//		double trendBeg = trend[trend.length - periodicity - 1];
//		double trendSlope = trendEnd - trendBeg;
//		trendSlope /= periodicity;

		for (int i = 0; i < periodicity; ++i) {
			int fromIndex = seasonal.length - periodicity + i; // [length - periodicity, length - 1]
			modelSeasonal[i] = seasonal[fromIndex];
			int fromExtended = seasonal.length + periodicity + i; // [seasonal.length + periodicity, seasonal.length + 2 * periodicity - 1]
			modelSeasonal[i + periodicity] = extendedSeasonal[fromExtended];
//			modelTrend[i] = trend[fromIndex];
//			modelTrend[i + periodicity] = trendEnd + (i + 1) * trendSlope;
		}

		LoessSettings settings = new LoessSettings(13, 2, 1);
		LoessSmoother smoother = new LoessSmoother(settings, modelSeasonal, null);

		double[] smoothedModelSeasonal = smoother.smooth();

		SummaryStatistics baseStats = new SummaryStatistics();
		SummaryStatistics smoothedStats = new SummaryStatistics();
		for (int i = 0; i < modelSeasonal.length - 1; ++i) {
			baseStats.addValue(modelSeasonal[i + 1] - modelSeasonal[i]);
			smoothedStats.addValue(smoothedModelSeasonal[i + 1] - smoothedModelSeasonal[i]);
		}

		assertTrue(baseStats.getMin() < smoothedStats.getMin());
		assertTrue(baseStats.getMax() > smoothedStats.getMax());
		assertTrue(baseStats.getStandardDeviation() > 0.5 * smoothedStats.getStandardDeviation());
	}

	@Ignore("Not really a test - used to generate data for comparison with python")
	@Test
	public void generateSmoothedDataTest() throws IOException {
		final double[] data = fTestData.values;

		int periodicity = 168;
		SeasonalTrendLoess stl = new SeasonalTrendLoess(data, periodicity, 2001, 1, 15);

		stl.decompose();

		TestUtilities.dumpStlResultsToFile(data, stl, "/tmp/stl_java.csv");

		stl.smoothSeasonal(13);

		TestUtilities.dumpStlResultsToFile(data, stl, "/tmp/stl_java_smoothed.csv");
	}
}

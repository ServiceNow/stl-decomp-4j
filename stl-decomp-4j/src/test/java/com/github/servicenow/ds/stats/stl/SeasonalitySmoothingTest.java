package com.github.servicenow.ds.stats.stl;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Ignore;
import org.junit.Test;

public class SeasonalitySmoothingTest {

	public static final double EPS = 1.0e-15;
	private final StlTestDataGenerator fTestData = new StlTestDataGenerator();

	@Test
	public void smoothSeasonalityTest() throws IOException {

		final double[] data = fTestData.values;

		int periodicity = 168;

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(periodicity).setSeasonalWidth(2001);
		builder.setInnerIterations(1).setRobustnessIterations(15);

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		double epsilon = 5.0e-03; // Doing absolute differences and some values > 10

		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResidual();

		for (int i = 0; i < data.length; ++i) {
			assertEquals(String.format("seasonal[%d]", i), fTestData.trend[i], trend[i], epsilon);
			assertEquals(String.format("trend[%d]", i), fTestData.seasonal[i], seasonal[i], 20 * epsilon);
			assertEquals(String.format("residuals[%d]", i), fTestData.residuals[i], residuals[i],
					20 * epsilon);
		}

		double[] extendedSeasonal = new double[seasonal.length + 2 * periodicity];

		final CyclicSubSeriesSmoother.Builder csBuilder = new CyclicSubSeriesSmoother.Builder().setWidth(2001);
		csBuilder.setDataLength(seasonal.length).setPeriodicity(periodicity).extrapolateForwardAndBack(1);
		CyclicSubSeriesSmoother seasonal_smoother = csBuilder.build();

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

		final LoessSmoother.Builder loessBuilder = new LoessSmoother.Builder().setWidth(13).setDegree(2).setJump(1);
		LoessSmoother smoother = loessBuilder.setData(modelSeasonal).build();

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

	@Test
	public void seasonalSmootherMinimalWidthTest() {
		final double[] data = fTestData.values;

		int periodicity = 168;

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(periodicity).setSeasonalWidth(2001);
		builder.setInnerIterations(1).setRobustnessIterations(15);

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		double[] trend = stl.getTrend().clone();
		double[] seasonal = stl.getSeasonal().clone(); // Make a copy
		double[] residuals = stl.getResidual().clone();

		stl.smoothSeasonal(3, false);

		for (int i = 1; i < seasonal.length - 1; ++i) {
			assertEquals("Smoothing with width 3 should have no effect", seasonal[i], stl.getSeasonal()[i], EPS);
			assertEquals("Smoothing with width 3 should have no effect", trend[i], stl.getTrend()[i], EPS);
			assertEquals("Smoothing with width 3 should have no effect", residuals[i], stl.getResidual()[i], 10 * EPS);
		}

		assertNotEquals(seasonal[0], stl.getSeasonal()[0], EPS);
		assertNotEquals(seasonal[seasonal.length - 1], stl.getSeasonal()[seasonal.length - 1], EPS);
	}

	@Test
	public void seasonalSmootherMinimalWidthNoEndpointFixTest() {
		final double[] data = fTestData.values;

		int periodicity = 168;

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(periodicity).setSeasonalWidth(2001);
		builder.setInnerIterations(1).setRobustnessIterations(15);

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		double[] trend = stl.getTrend().clone();
		double[] seasonal = stl.getSeasonal().clone(); // Make a copy
		double[] residuals = stl.getResidual().clone();

		stl.smoothSeasonal(3, true);

		assertArrayEquals("Smoothing with width 3 should have no effect", seasonal, stl.getSeasonal(), EPS);
		assertArrayEquals("Smoothing with width 3 should have no effect", trend, stl.getTrend(), EPS);
		assertArrayEquals("Smoothing with width 3 should have no effect", residuals, stl.getResidual(), 10 * EPS);
	}

	@Test
	public void seasonalSmootherWidth4Test() {
		SeasonalTrendLoess.Decomposition stl5 = getTestDecompositionWithSmootherWidth(5);
		SeasonalTrendLoess.Decomposition stl4 = getTestDecompositionWithSmootherWidth(4);

		compareDecompositions("Width 4 should be reset to 5", stl5, stl4);
	}

	@Test
	public void seasonalSmootherWidth2Test() {
		SeasonalTrendLoess.Decomposition stl3 = getTestDecompositionWithSmootherWidth(3);
		SeasonalTrendLoess.Decomposition stl2 = getTestDecompositionWithSmootherWidth(2);

		compareDecompositions("Width 2 should be reset to 3", stl3, stl2);
	}

	@Test
	public void seasonalSmootherWidth1Test() {
		SeasonalTrendLoess.Decomposition stl3 = getTestDecompositionWithSmootherWidth(3);
		SeasonalTrendLoess.Decomposition stl1 = getTestDecompositionWithSmootherWidth(1);

		compareDecompositions("Width 1 should be reset to 3", stl3, stl1);
	}

	@Test
	public void seasonalSmootherWidth0Test() {
		SeasonalTrendLoess.Decomposition stl3 = getTestDecompositionWithSmootherWidth(3);
		SeasonalTrendLoess.Decomposition stl0 = getTestDecompositionWithSmootherWidth(0);

		compareDecompositions("Width 0 should be reset to 3", stl3, stl0);
	}

	private void compareDecompositions(
			String message, SeasonalTrendLoess.Decomposition expected, SeasonalTrendLoess.Decomposition actual) {
		assertArrayEquals(message, expected.getSeasonal(), actual.getSeasonal(), EPS);
		assertArrayEquals(message, expected.getTrend(), actual.getTrend(), EPS);
		assertArrayEquals(message, expected.getResidual(), actual.getResidual(), 10 * EPS);
	}

	private SeasonalTrendLoess.Decomposition getTestDecompositionWithSmootherWidth(int width) {
		final double[] data = fTestData.values;

		int periodicity = 168;

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(periodicity).setSeasonalWidth(2001);
		builder.setInnerIterations(1).setRobustnessIterations(15);

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		stl.smoothSeasonal(width, true);
		return stl;
	}

	@Ignore("Not really a test - used to generate data for comparison with python")
	@Test
	public void generateSmoothedDataTest() throws IOException {
		final double[] data = fTestData.values;

		int periodicity = 168;

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(periodicity).setSeasonalWidth(2001);
		builder.setInnerIterations(1).setRobustnessIterations(15);

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		TestUtilities.dumpStlResultsToFile(data, stl, "/tmp/stl_java.csv");

		stl.smoothSeasonal(13);

		TestUtilities.dumpStlResultsToFile(data, stl, "/tmp/stl_java_smoothed.csv");
	}
}

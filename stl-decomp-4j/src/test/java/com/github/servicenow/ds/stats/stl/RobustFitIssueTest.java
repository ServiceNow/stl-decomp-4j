package com.github.servicenow.ds.stats.stl;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class RobustFitIssueTest {

	@Test
	public void testTwoWeekPeriodicFitNonRobust() {
		// With only two periods of data, I expect the periodic non-robust result to just give
		// the seasonal + trend equal to the average of the two sub-cycle values that are used
		// to compute the seasonal component. This test verifies that this is true.

		double[] data = SimulatedWeeklyMetric.getValues();

		assert data.length == 2016;

		int periodicity = 1008; // Exactly half the data

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();

		builder.setPeriodLength(periodicity)
				.setPeriodic()
				.setFlatTrend();

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		double[] seasonal = stl.getSeasonal().clone();

		assertPeriodic("Seasonal should be periodic", seasonal, periodicity, periodicity, 1.0e-15);

		double base = stl.getTrend()[0];

		for (int i = 0; i < periodicity; ++i) {
			double expected = (data[i] + data[i + periodicity]) / 2.0;
			double modeled = base + seasonal[i];
			assertEquals(expected, modeled, 5.0e-14);
		}
	}

	@Test
	public void testTwoWeekLinearFitNonRobust() {
		// With only two periods of data, the non-robust result with linear variability
		// in the seasonality and trend to overfit the data, resulting in zero residual.

		// Note that this will not happen with flat trend as the seasonality has to have
		// zero average over the entire two periods and the data won't necessarily, so
		// a linear trend is required to allow this and completely overfit the input.

		double[] data = SimulatedWeeklyMetric.getValues();

		assert data.length == 2016;

		int periodicity = 1008; // Exactly half the data

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();

		builder.setPeriodLength(periodicity)
				.setSeasonalWidth(100 * data.length)
				.setSeasonalDegree(1);

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		for (int i = 0; i < data.length; ++i)
			assertEquals(0.0, stl.getResidual()[i], 1.0e-13);
	}

	@Test
	public void testFourWeekRobustFitOutliers() {
		double[] data = SimulatedWeeklyMetric.getFourWeekValues();

		assert data.length == 4032;

		int periodicity = 1008;

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();

		builder.setPeriodLength(periodicity)
				.setRobust()
				.setPeriodic()
				.setFlatTrend();

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();


		double[] weights = stl.getWeights();

		int allZeroWeightCount = 0;
		for (int i = 0; i < periodicity; ++i) {
			boolean all = true;
			for (int j = 0; j < 4; ++j)
				all = all && weights[i + j * periodicity] == 0.0;

			if (all)
				++allZeroWeightCount;
		}

		assertEquals(0, allZeroWeightCount);

//		stl.smoothSeasonal(73);
//		printStlFit(stl);
	}

	@Test
	public void testTwoWeekRobustFitOutliers() {
		// This test demonstrates a problem with using robust fitting with only two periods.
		// The test is set up to pass but the upshot is that this setup is dangerous and should not be used.

		double[] data = SimulatedWeeklyMetric.getValues();

		assert data.length == 2016;

		int periodicity = 1008; // Exactly half the data

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();

		builder.setPeriodLength(periodicity)
				.setRobust()
				.setPeriodic()
				.setFlatTrend();

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		double[] seasonal = stl.getSeasonal().clone();
		double trend = stl.getTrend()[0];
		double[] weights = stl.getWeights();

		// There are two seasonal sub-series that have all the weights == 0, which basically blows
		// up the robustness iteration. When the number of iterations is odd, we find that the
		// weights come out zero and the points are set to the data points, which of course are not
		// periodic, so periodicity is badly broken.

		assertEquals(0.0, weights[497], 1.0e-15);
		assertEquals(0.0, weights[497 + periodicity], 1.0e-15);
		assertEquals(data[497], seasonal[497] + trend, 1.0e-2);
		assertEquals(data[497 + periodicity], seasonal[497 + periodicity] + trend, 1.0e-2);

		assertEquals(0.0, weights[927], 1.0e-15);
		assertEquals(0.0, weights[927 + periodicity], 1.0e-15);
		assertEquals(data[927], seasonal[927] + trend, 1.0e-2);
		assertEquals(data[927 + periodicity], seasonal[927 + periodicity] + trend, 1.0e-2);

		// Check that these are the only two places for which both weights are zero

		int zeros = 0;
		for (int i = 0; i < weights.length / 2; ++i)
			if (weights[i] == 0.0 && weights[i + periodicity] == 0.0)
				++zeros;

		assertEquals(2, zeros);

		// Aside from these two points, the seasonality is approximately periodic:

		seasonal[497] = seasonal[497 + periodicity];
		seasonal[927] = seasonal[927 + periodicity];
		assertPeriodic("Seasonal should be periodic", seasonal, periodicity, periodicity, 1.0e-2);

		// If we do an even number of iterations, then the situation is completely different - instead the solution
		// jumps back to something close to the non-robust solution.

		builder.setRobustnessIterations(16); // One more iteration than default.

		stlSmoother = builder.buildSmoother(data);

		stl = stlSmoother.decompose();

		seasonal = stl.getSeasonal().clone();
		trend = stl.getTrend()[0];
		weights = stl.getWeights();

		assertEquals(1.0, weights[497], 1.0e-4);
		assertEquals(1.0, weights[497 + periodicity], 1.0e-4);
		assertEquals(1.0, weights[927], 1.0e-4);
		assertEquals(1.0, weights[927 + periodicity], 1.0e-4);

		assertEquals((data[497] + data[497 + periodicity]) / 2, seasonal[497] + trend, 1.0e-2);
		assertEquals((data[927] + data[927 + periodicity]) / 2, seasonal[927] + trend, 1.0e-2);
	}

//	@Test
	public void printResidualConvergence() {
		double[] data = SimulatedWeeklyMetric.getValues();

		assert data.length == 2016;

		int periodicity = 1008; // Exactly half the data

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();

		builder.setPeriodLength(periodicity)
				.setRobust()
				.setPeriodic()
				.setFlatTrend();

		for (int i = 0; i < 100; ++i) {
			printForRobustnessIterations(i, data, builder);
		}
	}

	private static void assertPeriodic(String message, double[] array, int offset, int length, double delta) {
		for (int i = 0; i < length; ++i) {
			double expected = array[i];
			double actual = array[i + offset];

			if (doubleIsDifferent(expected, actual, delta))
				failNotEquals(message, i, i + offset, expected, actual);
		}
	}

	private static boolean doubleIsDifferent(double d1, double d2, double delta) {
		if (Double.compare(d1, d2) == 0) {
			return false;
		} else {
			return Math.abs(d1 - d2) > delta;
		}
	}

	private static void failNotEquals(String message, int first, int second, double expected, double actual) {
		String formatted = String.format(
				"%s: Expected[%d]: %f, Actual[%d] %f", message, first, expected, second, actual);
		throw new AssertionError(formatted);
	}

	private void printForRobustnessIterations(int iterations, double[] data, SeasonalTrendLoess.Builder builder) {

//		System.out.println(String.format("Analyzing for %d robustness iterations", iterations));

		builder.setRobustnessIterations(iterations);

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

//		printStlResultsAtIndex(497, data, stl);
//		printStlResultsAtIndex(497 + 1008, data, stl);
//		printStlResultsAtIndex(927, data, stl);
//		printStlResultsAtIndex(927 + 1008, data, stl);
//		printStlResultsAtIndex(2, data, stl);
//		printStlResultsAtIndex(1010, data, stl);

		printMeanSquareResidual(data, stl);
	}

	private void printStlResultsAtIndex(int idx, double[] data, SeasonalTrendLoess.Decomposition stl) {
		System.out.println(String.format("@ %4d: value = %f, trend = %f, seasonal = %f, residual = %f, weight = %f",
				idx, data[idx], stl.getTrend()[idx], stl.getSeasonal()[idx], stl.getResidual()[idx], stl.getWeights()[idx]));
	}

	private void printMeanSquareResidual(double[] data, SeasonalTrendLoess.Decomposition stl) {
		double r2 = 0;
		double w2 = 0;
		double[] residual = stl.getResidual();
		double[] weights = stl.getWeights();
		double r2max = 0;
		double wmin = 1;
		int idx_wmin = 0;
		double wmax = 0;
		int idx_wmax = 0;
		for (int i = 0; i < residual.length; ++i) {
			double res = residual[i];
			r2 += res * res;
			double w = 1 - weights[i];
			w2 += w * w;

			if (r2 > r2max)
				r2max = r2;

			if ((1 - w) < wmin) {
				wmin = (1 - w);
				idx_wmin = i;
			}

			if ((1 - w) > wmax && (1 - w) < .8) {
				wmax = 1 - w;
				idx_wmax = i;
			}
		}

		System.out.println(String.format("%f, %f, %f, %f, %d, %f , %f, %d, %f",
				Math.sqrt(r2 / data.length), r2max, Math.sqrt(w2 / data.length),
				wmin, idx_wmin, weights[376], wmax, idx_wmax, weights[816]));
	}

	private void printStlFit(SeasonalTrendLoess.Decomposition stl) {
		double[] seasonal = stl.getSeasonal();
		double[] trend = stl.getTrend();
		double[] residual = stl.getResidual();
		double[] weight = stl.getWeights();

		System.out.println("seasonal, trend, residual, weight");
		for (int i = 0; i < seasonal.length; ++i)
			System.out.println(String.format("%f, %f, %f, %f", seasonal[i], trend[i], residual[i], weight[i]));

	}
}

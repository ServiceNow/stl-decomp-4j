package com.github.servicenow.ds.stats.stl;

import org.junit.Ignore;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

/**
 * Integration tests for SeasonalTrendLoess (STL) decomposition.
 *
 * Created by Jim Crotinger on 20-Apr-2016.
 */

public class SeasonalTrendLoessTest {

	private final StlTestDataGenerator testDataGenerator = new StlTestDataGenerator();

	@Test
	public void pureSineTest() {

		long seed = 1234567L; // System.nanoTime() // change this to do random stress test

		double[] data = testDataGenerator.createNoisySeasonalData(144, 12, 1.0, 0.0, 0.0, seed);

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12).setSeasonalWidth(7).setNonRobust();
		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResidual();

		for (int i = 0; i < data.length; ++i) {
			assertEquals(data[i], seasonal[i], 1.0e-14);
			assertEquals(0.0, trend[i], 1.0e-14);
			assertEquals(0.0, residuals[i], 1.0e-14);
		}
	}

	@Test
	public void pureTrendTest() {

		long seed = 1234567L; // System.nanoTime() // change this to do random stress test

		double[] data = testDataGenerator.createNoisySeasonalData(144, 12, 0.0, 1.0, 0.0, seed);

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12).setSeasonalWidth(7).setNonRobust();
		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResidual();
		double[] weights = stl.getWeights();

		for (int i = 0; i < data.length; ++i) {
			assertEquals(data[i], trend[i], 1.0e-12);
			assertEquals(0.0, seasonal[i], 1.0e-12);
			assertEquals(0.0, residuals[i], 1.0e-12);
			assertEquals(1.0, weights[i], 1.0e-13);
		}
	}

	@Test
	public void squareWaveTest() {

		double[] data = testDataGenerator.createSquareWaveData();

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(288).setSeasonalWidth(13).setNonRobust();
		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResidual();
		double[] weights = stl.getWeights();

		for (int i = 0; i < data.length; ++i) {
			assertEquals(42.5, trend[i], 1.0e-12);
			assertEquals(data[i], seasonal[i] + stl.getTrend()[i], 1.0e-12);
			assertEquals(0.0, residuals[i], 1.0e-12);
			assertEquals(1.0, weights[i], 1.0e-13);
		}
	}

	@Test
	public void nonRobustRegressionTest() {

		// Check results against a Python/Fortran test on noisy/trendy/sinusoidal data.

		double[] data = new double[144];
		for (int i = 0; i < data.length; ++i) {
			data[i] = fNonRobustNoisySinusoidResults[i][0];
		}

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12).setSeasonalWidth(7);
		builder.setInnerIterations(2).setRobustnessIterations(0); // Should be the same as setNonRobust()
		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		//dumpStlResultsToFile(data, stl, "/tmp/stl_non_robust_test.csv");

		double epsilon = 7.0e-11; // Doing absolute differences and some values > 10

		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResidual();
		double[] weights = stl.getWeights();

		for (int i = 0; i < data.length; ++i) {
			assertEquals(String.format("trend[%d]", i), fNonRobustNoisySinusoidResults[i][1], trend[i], epsilon);
			assertEquals(String.format("seasonal[%d]", i), fNonRobustNoisySinusoidResults[i][2], seasonal[i], epsilon);
			assertEquals(String.format("residuals[%d]", i), fNonRobustNoisySinusoidResults[i][3], residuals[i], epsilon);
			assertEquals(1.0, weights[i], 1.0e-13);
		}
	}

	@Test
	public void forcedPeriodicityTest() {
		// To force periodicity, smooth the cyclic sub-series with a zero degree polynomial and a very long window,
		// causing it to interpolate each sub-series point at the average for the subseries.

		double[] data = new double[144];
		for (int i = 0; i < data.length; ++i) {
			data[i] = fNonRobustNoisySinusoidResults[i][0];
		}

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12);
		builder.setSeasonalWidth(100000001).setSeasonalDegree(0).setSeasonalJump(100000001); // Force periodic by hand
		builder.setTrendWidth(23);
		builder.setLowpassWidth(13);
		builder.setInnerIterations(2).setRobustnessIterations(0);

		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double epsilon = 0.0;

		double[] seasonal = stl.getSeasonal();

		for (int i = 0; i < 12; ++i) {
			for (int p = 0; p < 12; ++p) {
				assertEquals(seasonal[i], seasonal[i + p * 12], epsilon);
			}
		}
	}

	@Test
	public void setPeriodicTest() {
		// To force periodicity, smooth the cyclic sub-series with a zero degree polynomial and a very long window,
		// causing it to interpolate each sub-series point at the average for the subseries.

		double[] data = new double[144];
		for (int i = 0; i < data.length; ++i) {
			data[i] = fNonRobustNoisySinusoidResults[i][0];
		}

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12);
		builder.setPeriodic();
		builder.setTrendWidth(23);
		builder.setLowpassWidth(13);
		builder.setInnerIterations(2).setRobustnessIterations(0);

		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double epsilon = 2e-8; // setPeriodic only sets the length to 100 * data.length

		double[] seasonal = stl.getSeasonal();

		for (int i = 0; i < 12; ++i) {
			for (int p = 0; p < 12; ++p) {
				assertEquals(seasonal[i], seasonal[i + p * 12], epsilon);
			}
		}
	}

	@Test
	public void forcedPeriodicityTest2() {
		// Same as above but check results with different trend and lowpass settings.

		double[] data = new double[144];
		for (int i = 0; i < data.length; ++i) {
			data[i] = fNonRobustNoisySinusoidResults[i][0];
		}

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12);
		builder = builder.setSeasonalWidth(100000001).setSeasonalDegree(0).setSeasonalJump(100000001); // Force periodic by hand
		builder = builder.setTrendWidth(23).setTrendDegree(0).setTrendJump(100000);
		builder = builder.setLowpassWidth(13).setLowpassDegree(0).setLowpassJump(1);
		builder = builder.setInnerIterations(2).setRobustnessIterations(0);

		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double epsilon = 1.0e-10; // Doing absolute differences and some values > 10

		double[] seasonal = stl.getSeasonal();

		for (int i = 0; i < 12; ++i) {
			for (int p = 0; p < 12; ++p) {
				assertEquals(seasonal[i], seasonal[i + p * 12], epsilon);
			}
		}
	}

	@Test
	public void linearTrendTest() {
		long seed = 1234567L; // System.nanoTime() // change this to do random stress test

		double[] data = testDataGenerator.createNoisySeasonalData(144, 12, 1.0, 0.2, 0.1, seed);

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setLinearTrend().setPeriodLength(12).setSeasonalWidth(1000000).setRobust();

		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double[] trend = stl.getTrend();

		double dt = trend[1] - trend[0];
		for (int i = 1; i < trend.length; ++i)
			assertEquals(dt, trend[i] - trend[i-1], 1.0e-14);

		double dx = 2 * Math.PI / 12;

		assertEquals(dt, 0.2 * dx, 1.0e-4);
	}


	@Test
	public void flatTrendTest() {
		long seed = 1234567L; // System.nanoTime() // change this to do random stress test

		double[] data = testDataGenerator.createNoisySeasonalData(144, 12, 1.0, 0.0, 0.1, seed);

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setFlatTrend().setPeriodLength(12).setSeasonalWidth(1000000).setRobust();

		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double[] trend = stl.getTrend();

		for (int i = 1; i < trend.length; ++i)
			assertEquals(0.0, trend[i] - trend[i-1], 1.0e-13);
	}

	@Test
	public void sineWithOutlierTest() {

		long seed = 1234567L; // System.nanoTime() // change this to do random stress test

		double[] data = testDataGenerator.createNoisySeasonalData(144, 12, 1.0, 0.0, 0.0, seed);
		data[100] = 1000;

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12).setSeasonalWidth(1000000).setRobust();

		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		// printStlResults(data, stl);

		double epsilon = 1.0e-4;

		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResidual();

		for (int i = 0; i < data.length; ++i) {
			if (i != 100) {
				assertEquals(String.format("seasonal[%d]", i), data[i], seasonal[i], epsilon);
				assertEquals(String.format("trend[%d]", i), 0.0, trend[i], epsilon);
				assertEquals(String.format("residuals[%d]", i), 0.0, residuals[i], epsilon);
			} else {
				assertEquals(String.format("seasonal[%d]", i), data[i - 12], seasonal[i], epsilon);
				assertEquals(String.format("trend[%d]", i), 0.0, trend[i], epsilon);
				assertEquals(String.format("residuals[%d]", i), 1.0, residuals[i] / 1000.0, 1.0e-3);
			}
		}
	}

	@Test
	public void robustRegressionTest() {

		// Check results against a Python/Fortran test on noisy/trendy/sinusoidal data.

		double[] data = new double[144];
		for (int i = 0; i < data.length; ++i) {
			data[i] = fRobustNoisySinusoidResults[i][0];
		}

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12).setSeasonalWidth(7);
		builder.setInnerIterations(1).setRobustnessIterations(1);

		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		// TODO: The differences here are larger than I had expected. They're quite small but significantly larger than
		// machine epsilon, and looking at the difference between the seasonality from Java and Python/Fortran shows
		// some spurious stuff at the right boundary. Investigate in the future.

		double epsilon = 2.0e-07; // Doing absolute differences and some values > 10

		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResidual();

		for (int i = 0; i < data.length; ++i) {
			assertEquals(String.format("seasonal[%d]", i), fRobustNoisySinusoidResults[i][1], trend[i], epsilon);
			assertEquals(String.format("trend[%d]", i), fRobustNoisySinusoidResults[i][2], seasonal[i], epsilon);
			assertEquals(String.format("residuals[%d]", i), fRobustNoisySinusoidResults[i][3], residuals[i],
					epsilon);
		}
	}

	@Test
	public void periodicBuilderCanBeReused() {
		double[] data = SimulatedWeeklyMetric.getFourWeekValues();

		int periodicity = 1008;

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();

		builder.setPeriodLength(periodicity)
				.setRobust()
				.setPeriodic()
				.setFlatTrend();

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		assertNotNull(stl);

		builder.setRobustnessIterations(17);

		// Previously this resulted in an exception:
		// SeasonalTrendLoess.Builder: setSeasonalWidth and setPeriodic cannot both be called.

		stlSmoother = builder.buildSmoother(data);

		stl = stlSmoother.decompose();

		assertNotNull(stl);
	}

	@Test
	public void linearTrendBuilderCanBeReused() {
		double[] data = SimulatedWeeklyMetric.getFourWeekValues();

		int periodicity = 1008;

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();

		builder.setPeriodLength(periodicity)
				.setRobust()
				.setLinearTrend()
				.setSeasonalWidth(101);

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		assertNotNull(stl);

		builder.setRobustnessIterations(17);

		// Previously this resulted in an exception:
		// SeasonalTrendLoess.Builder: setTrendWidth incompatible with flat/linear trend.

		stlSmoother = builder.buildSmoother(data);

		stl = stlSmoother.decompose();

		assertNotNull(stl);
	}

	@Test
	public void flatTrendBuilderCanBeReused() {
		double[] data = SimulatedWeeklyMetric.getFourWeekValues();

		int periodicity = 1008;

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();

		builder.setPeriodLength(periodicity)
				.setRobust()
				.setFlatTrend()
				.setSeasonalWidth(101);

		SeasonalTrendLoess stlSmoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = stlSmoother.decompose();

		assertNotNull(stl);

		builder.setRobustnessIterations(17);

		// Previously this resulted in an exception:
		// SeasonalTrendLoess.Builder: setTrendWidth incompatible with flat/linear trend.

		stlSmoother = builder.buildSmoother(data);

		stl = stlSmoother.decompose();

		assertNotNull(stl);
	}

	@Test(expected=IllegalArgumentException.class)
	public void periodicityMustBeAtLeastTwo() {
		new SeasonalTrendLoess.Builder().setPeriodLength(1);
	}

	@Test(expected=IllegalArgumentException.class)
	public void dataMustHaveAtLeastTwoPeriods() {
		double[] data = testDataGenerator.createNoisySeasonalData(144, 12, 1.0, 0.0, 0.0, 123L);
		new SeasonalTrendLoess.Builder().setPeriodLength(120).setSeasonalWidth(999).setNonRobust().buildSmoother(data);
	}

	@Test(expected=IllegalArgumentException.class)
	public void nullDataThrows() {
		new SeasonalTrendLoess.Builder().setPeriodLength(120).setSeasonalWidth(999).buildSmoother(null);
	}

	@Test(expected=IllegalArgumentException.class)
	public void seasonalWidthMustBeSet() {
		new SeasonalTrendLoess.Builder().setPeriodLength(120).buildSmoother(new double[2000]);
	}

	@Test(expected=IllegalArgumentException.class)
	public void periodLengthMustBeSet() {
		new SeasonalTrendLoess.Builder().setSeasonalWidth(999).buildSmoother(new double[2000]);
	}

	@Test(expected=IllegalArgumentException.class)
	public void setPeriodicDisallowsSeasonalWidth() {
		getTestBuilder().setSeasonalWidth(999).buildSmoother(new double[2000]);
	}

	@Test(expected=IllegalArgumentException.class)
	public void setPeriodicDisallowsSeasonalJump() {
		getTestBuilder().setSeasonalDegree(2).buildSmoother(new double[2000]);
	}

	@Test(expected=IllegalArgumentException.class)
	public void setPeriodicDisallowsSeasonalDegree() {
		getTestBuilder().setSeasonalJump(1).buildSmoother(new double[2000]);
	}

	@Test(expected=IllegalArgumentException.class)
	public void setFlatTrendDisallowsTrendWidth() {
		getTestBuilder().setFlatTrend().setTrendWidth(999).buildSmoother(new double[2000]);
	}

	@Test(expected=IllegalArgumentException.class)
	public void setFlatTrendDisallowsTrendJump() {
		getTestBuilder().setFlatTrend().setTrendJump(1).buildSmoother(new double[2000]);
	}

	@Test(expected=IllegalArgumentException.class)
	public void setFlatTrendDisallowsTrendDegree() {
		getTestBuilder().setFlatTrend().setTrendDegree(2).buildSmoother(new double[2000]);
	}

	@Test(expected=IllegalArgumentException.class)
	public void setLinearTrendDisallowsTrendWidth() {
		getTestBuilder().setLinearTrend().setTrendWidth(999).buildSmoother(new double[2000]);
	}

	@Test(expected=IllegalArgumentException.class)
	public void setLinearTrendDisallowsTrendJump() {
		getTestBuilder().setLinearTrend().setTrendJump(1).buildSmoother(new double[2000]);
	}

	@Test(expected=IllegalArgumentException.class)
	public void setLinearTrendDisallowsTrendDegree() {
		getTestBuilder().setLinearTrend().setTrendDegree(2).buildSmoother(new double[2000]);
	}

	private SeasonalTrendLoess.Builder getTestBuilder() {
		return new SeasonalTrendLoess.Builder().setPeriodic().setPeriodLength(10);
	}

	@Ignore("For manual testing only")
	@Test
	public void pureNoiseTest() {

		double[] data = testDataGenerator.createNoisySeasonalDataWithTimeSeed(144, 12, 0.0, 0.0, 1.0);

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12).setSeasonalWidth(7).setRobust();

		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResidual();

		double rsum = 0.0;
		double r2sum = 0.0;
		double dsum = 0.0;
		double d2sum = 0.0;
		double tsum = 0.0;
		double t2sum = 0.0;
		double ssum = 0.0;
		double s2sum = 0.0;
		for (int i = 0; i < data.length; ++i) {
			final double r = residuals[i];
			rsum += r;
			r2sum += r * r;

			final double d = data[i];
			dsum += d;
			d2sum += d * d;

			final double t = trend[i];
			tsum += t;
			t2sum += t * t;

			final double s = seasonal[i];
			ssum += s;
			s2sum += s * s;
		}

		double trendMean = tsum / data.length;
		double trendVar = (t2sum - trendMean * trendMean) / (data.length - 1);

		double resMean = rsum / data.length;
		double resVar = (r2sum - resMean * resMean) / (data.length - 1);

		double dataMean = dsum / data.length;
		double dataVar = (d2sum - dataMean * dataMean) / (data.length - 1);

		double seasonalMean = ssum / data.length;
		double seasonalVar = (s2sum - seasonalMean * seasonalMean) / (data.length - 1);

		final double stlMean = seasonalMean + trendMean + resMean;
		assertEquals(dataMean, stlMean, 1.0e-8);
		final double stlVar = seasonalVar + trendVar + resVar;

		System.out.println("     \t\tmean\tvariance");
		System.out.printf("data     \t%f\t%f%n", dataMean, dataVar);
		System.out.printf("seasonal \t%f\t%f%n", seasonalMean, seasonalVar);
		System.out.printf("trend    \t%f\t%f%n", trendMean, trendVar);
		System.out.printf("residual \t%f\t%f%n", resMean, resVar);
		System.out.printf("stl sum  \t%f\t%f%n", stlMean, stlVar);
	}

	@Ignore("This is just for data generation so we can compare the Java results in Python")
	@Test
	public void stlDataFactory() {

		double[] data = testDataGenerator.createNoisySeasonalDataWithTimeSeed(144, 12, 10.0, 1.0, 2.0);

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12).setSeasonalWidth(7).setNonRobust();

		SeasonalTrendLoess smoother = builder.buildSmoother(data);

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		double[] residuals = stl.getResidual();

		double rsum = 0.0;
		double r2sum = 0.0;
		for (int i = 0; i < data.length; ++i) {
			double r = residuals[i];
			rsum += r;
			r2sum += r * r;
		}

		double mean = rsum / data.length;
		double variance = (r2sum - rsum * rsum) / (data.length - 1);

		System.out.printf("Residual has mean %f and variance %f%n", mean, variance);

		// dumpStlResultsToFile(data, stl, "/tmp/stl_test.csv");
	}

	@Test
	public void toStringTest() {
		double[] data = testDataGenerator.createNoisySeasonalDataWithTimeSeed(144, 12, 10.0, 1.0, 2.0);

		SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder().setPeriodLength(12).setSeasonalWidth(7).setNonRobust();

		SeasonalTrendLoess stl = builder.buildSmoother(data);

		assertEquals(
				"SeasonalTrendLoess: [\n" +
				"inner iterations     = 2\n" +
				"outer iterations     = 0\n" +
				"periodicity          = 12\n" +
				"seasonality settings = [width = 7, degree = 1, jump = 1]\n" +
				"trend settings       = [width = 23, degree = 1, jump = 3]\n" +
				"lowpass settings     = [width = 13, degree = 1, jump = 2]\n]", stl.toString());
	}

	@SuppressWarnings("unused")
	private void printStlResults(double[] data, SeasonalTrendLoess.Decomposition stl) {
		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResidual();
		double[] weights = stl.getWeights();

		System.out.printf("%8s %18s  \t%18s  \t%18s  \t%18s  \t%18s%n", "index", "data", "trend",
				"seasonal", "residual", "weights");
		for (int i = 0; i < data.length; ++i) {
			System.out.printf("%8d %18.14f  \t%18.14f  \t%18.14f  \t%18.14f  \t%18.14f%n", i, data[i],
					trend[i], seasonal[i], residuals[i], weights[i]);
		}
	}

	/**
	 * The arrays below contain results from running STL in Python.
	 */
	private final double[][] fNonRobustNoisySinusoidResults = new double[][] {
			{ 1.34006384538, 1.05868325038, -0.637510654892, 0.918891249895, 1.0 },
			{ 3.97126146325, 1.36552711582, 2.12635406017, 0.479380287255, 1.0 },
			{ 11.9840016758, 1.67237098126, 11.2169231848, -0.905292490221, 1.0 },
			{ 13.9575483181, 1.97921484671, 11.5329348655, 0.445398605914, 1.0 },
			{ 13.5012684434, 2.31316292988, 9.16795160766, 2.0201539059, 1.0 },
			{ 5.92000248407, 2.64711101304, 4.0621611416, -0.789269670575, 1.0 },
			{ 4.04443595375, 2.98105909621, 1.28365525234, -0.220278394806, 1.0 },
			{ -2.31984699321, 3.33913859485, -5.42126778914, -0.237717798913, 1.0 },
			{ -7.44652983694, 3.69721809349, -9.87532493022, -1.2684230002, 1.0 },
			{ -7.59922481848, 4.05529759212, -9.62466694595, -2.02985546465, 1.0 },
			{ -5.81698041276, 4.4829095444, -9.06664360027, -1.23324635689, 1.0 },
			{ 1.18349327913, 4.91052149667, -4.23383090777, 0.50680269023, 1.0 },
			{ 3.42906032699, 5.33813344895, -0.216694005183, -1.69237911678, 1.0 },
			{ 8.86047912788, 5.95689889739, 2.29106917635, 0.612511054145, 1.0 },
			{ 16.6765050489, 6.57566434584, 10.5266672619, -0.425826558789, 1.0 },
			{ 17.266341472, 7.19442979428, 10.8855608097, -0.813649131972, 1.0 },
			{ 14.9436480866, 7.84012973324, 8.60686375448, -1.50334540117, 1.0 },
			{ 14.4078614747, 8.4858296722, 4.67168037241, 1.2503514301, 1.0 },
			{ 11.5767464114, 9.13152961117, 0.322283449786, 2.12293335044, 1.0 },
			{ 4.45438126235, 9.73410951694, -5.21657610629, -0.0631521483045, 1.0 },
			{ 1.95758451299, 10.3366894227, -9.17469848828, 0.795593578555, 1.0 },
			{ 2.81766795845, 10.9392693285, -9.28512677252, 1.16352540249, 1.0 },
			{ 3.43325948467, 11.477155394, -8.42519281898, 0.381296909681, 1.0 },
			{ 6.14152242692, 12.0150414594, -4.67504093985, -1.19847809268, 1.0 },
			{ 12.6224384274, 12.5529275249, 0.180622231164, -0.111111328719, 1.0 },
			{ 13.8305271556, 13.1386400654, 2.63954922775, -1.94766213759, 1.0 },
			{ 26.220792221, 13.7243526059, 9.75147688172, 2.74496273335, 1.0 },
			{ 24.9032272453, 14.3100651464, 10.284339572, 0.308822526888, 1.0 },
			{ 20.544524016, 14.9360998367, 8.1895380337, -2.58111385434, 1.0 },
			{ 20.0594121059, 15.562134527, 5.28710485031, -0.789827271383, 1.0 },
			{ 12.7303421157, 16.1881692173, -0.444541791158, -3.01328531044, 1.0 },
			{ 11.6043961988, 16.777540246, -5.10155202358, -0.0715920235923, 1.0 },
			{ 10.6824452466, 17.3669112747, -8.52631259072, 1.84184656271, 1.0 },
			{ 11.0727204154, 17.9562823033, -9.13056836665, 2.2470064787, 1.0 },
			{ 12.2983516401, 18.5277532344, -7.94467159281, 1.7152699986, 1.0 },
			{ 15.9523954008, 19.0992241654, -5.0486956052, 1.90186684063, 1.0 },
			{ 21.0970458051, 19.6706950964, 0.669502054744, 0.756848653873, 1.0 },
			{ 21.4027230908, 20.1602425619, 3.02026595341, -1.77778542454, 1.0 },
			{ 28.5814779548, 20.6497900275, 9.01623262609, -1.0845446987, 1.0 },
			{ 30.2996286764, 21.139337493, 9.82465660878, -0.664365425303, 1.0 },
			{ 30.3617189355, 21.5409897464, 7.98775330449, 0.832975884673, 1.0 },
			{ 28.1043671962, 21.9426419998, 5.75641778555, 0.405307410838, 1.0 },
			{ 20.4310624852, 22.3442942532, -1.30067379659, -0.612557971446, 1.0 },
			{ 18.4489246104, 22.8074562189, -5.15058056008, 0.79204895152, 1.0 },
			{ 13.8749491267, 23.2706181846, -7.98060752612, -1.41506153182, 1.0 },
			{ 14.7829297649, 23.7337801503, -9.25495754254, 0.30410715711, 1.0 },
			{ 17.0385310723, 24.2506542014, -7.67966572347, 0.467542594365, 1.0 },
			{ 16.9080305609, 24.7675282524, -5.20139814343, -2.65809954804, 1.0 },
			{ 26.9095584401, 25.2844023035, 0.771651107302, 0.853505029367, 1.0 },
			{ 32.805497603, 25.7930228806, 3.97208265814, 3.04039206427, 1.0 },
			{ 32.8213035832, 26.3016434578, 8.04987552411, -1.53021539867, 1.0 },
			{ 35.8117291252, 26.8102640349, 9.83741076587, -0.835945675643, 1.0 },
			{ 35.9188987252, 27.2873661796, 8.99865891469, -0.367126368998, 1.0 },
			{ 33.9233255585, 27.7644683242, 5.92117723818, 0.23767999615, 1.0 },
			{ 27.798615791, 28.2415704688, -0.971906950386, 0.52895227257, 1.0 },
			{ 24.6809488112, 28.7424106202, -6.03865696461, 1.97719515561, 1.0 },
			{ 21.3304330734, 29.2432507716, -7.94132333616, 0.0285056380376, 1.0 },
			{ 18.6458655751, 29.7440909229, -10.5423938971, -0.555831450704, 1.0 },
			{ 21.9038588092, 30.2333804574, -8.11192219276, -0.217599455381, 1.0 },
			{ 25.0923993495, 30.7226699918, -4.37345553337, -1.25681510895, 1.0 },
			{ 31.0479089569, 31.2119595263, 0.29306795158, -0.457118521014, 1.0 },
			{ 35.3072734235, 31.6665242388, 5.24357327285, -1.60282408818, 1.0 },
			{ 40.3046564168, 32.1210889514, 7.50762972686, 0.675937738626, 1.0 },
			{ 43.8903039089, 32.5756536639, 9.68149728408, 1.63315296089, 1.0 },
			{ 43.7848121552, 33.1486534991, 9.44913131297, 1.18702734308, 1.0 },
			{ 40.1001105393, 33.7216533344, 5.81595364981, 0.562503555057, 1.0 },
			{ 34.1394090546, 34.2946531696, -0.14649453212, -0.00874958296351, 1.0 },
			{ 24.7734384956, 34.9244170247, -6.61743800848, -3.53354052068, 1.0 },
			{ 28.4868264474, 35.5541808798, -7.85127188434, 0.78391745194, 1.0 },
			{ 23.5650392038, 36.1839447349, -11.0561083704, -1.5627971607, 1.0 },
			{ 26.1334345662, 36.7928083152, -8.30372678433, -2.35564696466, 1.0 },
			{ 36.401593842, 37.4016718954, -3.49124345754, 2.49116540414, 1.0 },
			{ 37.8406636071, 38.0105354757, -0.290376451724, 0.120504583078, 1.0 },
			{ 45.1336156778, 38.6977055658, 5.55950068875, 0.876409423266, 1.0 },
			{ 47.3840307948, 39.3848756559, 7.31413666921, 0.685018469625, 1.0 },
			{ 49.5183635921, 40.0720457461, 9.47800659334, -0.0316887473213, 1.0 },
			{ 50.4072881381, 40.7201825602, 9.24207092687, 0.445034650992, 1.0 },
			{ 46.1996378676, 41.3683193744, 5.68867773834, -0.857359245173, 1.0 },
			{ 42.2601415989, 42.0164561886, 0.295373264815, -0.0516878544983, 1.0 },
			{ 36.0938359413, 42.5136091039, -6.56403425207, 0.144261089425, 1.0 },
			{ 35.1415147161, 43.0107620192, -7.94273729353, 0.0734899903984, 1.0 },
			{ 32.0892994971, 43.5079149346, -10.7423375654, -0.676277872101, 1.0 },
			{ 38.7844760669, 43.9649645394, -8.37801127587, 3.19752280342, 1.0 },
			{ 42.7236303656, 44.4220141442, -3.15274839056, 1.45436461195, 1.0 },
			{ 43.1687334466, 44.879063749, -0.507720109852, -1.20261019252, 1.0 },
			{ 51.6229611522, 45.2917691913, 5.56681930685, 0.764372654079, 1.0 },
			{ 52.2252954846, 45.7044746336, 7.24699459784, -0.726173746882, 1.0 },
			{ 54.0910988294, 46.1171800759, 8.97641409926, -1.00249534578, 1.0 },
			{ 54.7234597154, 46.4443829234, 8.55742827258, -0.278351480577, 1.0 },
			{ 52.3844848676, 46.7715857708, 5.96360031469, -0.35070121789, 1.0 },
			{ 48.606331588, 47.0987886182, 0.341051445845, 1.16649152387, 1.0 },
			{ 41.7502176634, 47.4524274108, -5.24996709842, -0.452242648922, 1.0 },
			{ 39.5314465611, 47.8060662033, -8.81749687838, 0.54287723617, 1.0 },
			{ 40.59204504, 48.1597049958, -10.0662369284, 2.4985769726, 1.0 },
			{ 38.075991823, 48.5843162209, -8.47729779753, -2.03102660028, 1.0 },
			{ 43.6884728088, 49.0089274459, -3.53305933818, -1.78739529894, 1.0 },
			{ 50.2610716452, 49.433538671, -0.633269168473, 1.46080214271, 1.0 },
			{ 54.4960668534, 49.9028970783, 6.01917649429, -1.42600671919, 1.0 },
			{ 56.3803528362, 50.3722554856, 7.11868933079, -1.11059198015, 1.0 },
			{ 60.0838614108, 50.8416138929, 8.75332342601, 0.488924091919, 1.0 },
			{ 58.8956445746, 51.3863616045, 7.7520968821, -0.242813911977, 1.0 },
			{ 58.060931466, 51.9311093162, 6.73218381722, -0.602361667357, 1.0 },
			{ 51.7166992173, 52.4758570278, 0.340533050362, -1.09969086082, 1.0 },
			{ 50.5968330137, 53.0972270573, -4.40225818855, 1.90186414497, 1.0 },
			{ 44.1263683731, 53.7185970868, -9.64876105556, 0.0565323418103, 1.0 },
			{ 42.5421724941, 54.3399671163, -9.13838492738, -2.65940969482, 1.0 },
			{ 48.664593975, 54.9958206931, -9.47969949886, 3.14847278086, 1.0 },
			{ 50.8023147121, 55.6516742698, -4.25291185564, -0.596447702032, 1.0 },
			{ 54.5754039189, 56.3075278465, -0.687352796375, -1.0447711312, 1.0 },
			{ 62.4512633386, 56.92637349, 6.40880812044, -0.883918271832, 1.0 },
			{ 67.1464257363, 57.5452191335, 7.61306155329, 1.98814504959, 1.0 },
			{ 65.7645646262, 58.164064777, 9.10820192032, -1.50770207107, 1.0 },
			{ 65.8427013365, 58.7969504, 7.82898451133, -0.783233574781, 1.0 },
			{ 68.453524709, 59.429836023, 6.72172780171, 2.30196088431, 1.0 },
			{ 60.630284987, 60.062721646, 0.0673363516321, 0.500226989359, 1.0 },
			{ 57.5862119318, 60.6667634283, -4.4496856612, 1.36913416468, 1.0 },
			{ 48.7292236134, 61.2708052106, -9.15425914864, -3.3873224486, 1.0 },
			{ 55.1160500883, 61.8748469929, -9.11687821224, 2.35808130763, 1.0 },
			{ 47.9864061192, 62.3798706233, -9.82452238847, -4.56894211558, 1.0 },
			{ 60.6231791994, 62.8848942536, -5.39824766774, 3.13653261354, 1.0 },
			{ 63.6631990614, 63.389917884, -0.870914536818, 1.14419571427, 1.0 },
			{ 74.0142829683, 63.8762754731, 6.55657696159, 3.58143053358, 1.0 },
			{ 68.9249783922, 64.3626330623, 8.54109174975, -3.97874641984, 1.0 },
			{ 76.4283531358, 64.8489906515, 9.29041005625, 2.28895242811, 1.0 },
			{ 72.264141402, 65.3236142267, 8.08877390274, -1.1482467274, 1.0 },
			{ 73.1332843951, 65.7982378019, 6.57790552964, 0.757141063539, 1.0 },
			{ 66.2811787103, 66.2728613771, -0.212208779667, 0.220526112916, 1.0 },
			{ 58.8107417597, 66.7117400434, -4.29702400936, -3.60397427437, 1.0 },
			{ 59.4586934657, 67.1506187097, -8.63685308681, 0.944927842774, 1.0 },
			{ 58.248381867, 67.589497376, -9.07394410023, -0.267171408786, 1.0 },
			{ 60.5263995812, 68.0987003964, -10.4006322361, 2.82833142085, 1.0 },
			{ 61.1582405603, 68.6079034168, -6.48974270857, -0.959920148003, 1.0 },
			{ 67.4183160282, 69.1171064372, -1.07904324297, -0.61974716606, 1.0 },
			{ 74.4100207794, 69.6382259899, 6.78910117859, -2.0173063891, 1.0 },
			{ 81.9027381147, 70.1593455426, 9.45806244355, 2.28533012851, 1.0 },
			{ 79.118934638, 70.6804650953, 9.51965748011, -1.0811879374, 1.0 },
			{ 80.8959597864, 71.1877566465, 8.42727204701, 1.28093109288, 1.0 },
			{ 76.4955509714, 71.6950481976, 6.30604215845, -1.5055393847, 1.0 },
			{ 71.5154545294, 72.2023397487, -0.481831992404, -0.205053226907, 1.0 },
			{ 70.1564170333, 72.7105644776, -4.29770913288, 1.74356168852, 1.0 },
			{ 66.1529306932, 73.2187892065, -7.92211077793, 0.85625226465, 1.0 },
			{ 64.2350467473, 73.7270139354, -9.06147234997, -0.430494838099, 1.0 },
			{ 62.8779891919, 74.2340184888, -10.8613368773, -0.494692419587, 1.0 },
			{ 66.4852654904, 74.7410230421, -7.71329348619, -0.542464065513, 1.0 } };

	private final double[][] fRobustNoisySinusoidResults = new double[][] {
			{ -3.0189282813436, -0.0800823669536, -1.8687115842137, -1.0701343301763, 0.9412570824278 },
			{ 5.0778409942246, 0.5361271455863, 4.9127196692433, -0.3710058206049, 0.9967002109896 },
			{ 7.8291237106354, 1.1523366581261, 7.6955453238228, -1.0187582713135, 0.8717585393042 },
			{ 14.5176959771249, 1.7685461706660, 12.3074529116668, 0.4416968947921, 0.9857639268194 },
			{ 11.0097644574763, 2.3691976388519, 7.3903361388585, 1.2502306797659, 0.8864137929749 },
			{ 8.8754134904347, 2.9698491070377, 4.3550801429867, 1.5504842404103, 0.8530403314320 },
			{ 4.5347445838676, 3.5705005752236, 0.0715434568757, 0.8927005517683, 0.9460437270346 },
			{ -0.6070339602025, 4.1522309999183, -2.6045369599316, -2.1547280001891, 0.7604335360489 },
			{ -7.7208554791078, 4.7339614246129, -12.4324580845084, -0.0223588192123, 0.9944657737936 },
			{ -4.2241088734840, 5.3156918493076, -8.9428273943797, -0.5969733284119, 0.9637465871751 },
			{ 1.3947067273650, 5.8140026955752, -6.5277736154627, 2.1084776472526, 0.6856566596348 },
			{ 1.3266910351606, 6.3123135418428, -4.8567928525301, -0.1288296541521, 0.9902414899381 },
			{ 6.3219737880529, 6.8106243881105, -1.1685647458911, 0.6799141458335, 0.9941268843688 },
			{ 14.3218697206289, 7.2549644911939, 4.5383722267715, 2.5285330026634, 0.6974870916114 },
			{ 19.3976263092283, 7.6993045942774, 7.8544177119052, 3.8439040030457, 0.4747919526041 },
			{ 19.5139980433891, 8.1436446973609, 11.8638551331305, -0.4935017871022, 0.9825288000287 },
			{ 14.1508330850075, 8.5850467430852, 7.7610381784511, -2.1952518365288, 0.7439976272132 },
			{ 11.8046913258947, 9.0264487888095, 4.3217017060672, -1.5434591689820, 0.8606130972630 },
			{ 6.9127222031276, 9.4678508345338, -0.3363112470369, -2.2188173843694, 0.7599727322287 },
			{ 9.3797025107869, 9.8798316982575, -2.4766821693747, 1.9765529819041, 0.7801093396956 },
			{ -1.2987186194282, 10.2918125619811, -10.7299746406289, -0.8605565407805, 0.9369958308625 },
			{ 2.2419489466537, 10.7037934257048, -9.6532037514485, 1.1913592723974, 0.9204952837679 },
			{ 0.0141704379879, 11.2411047337361, -7.1932478872799, -4.0336864084683, 0.4045263986809 },
			{ 8.0447974841617, 11.7784160417674, -5.2770958705710, 1.5434773129653, 0.8707309361186 },
			{ 15.3540332761918, 12.3157273497987, -0.5383287655261, 3.5766346919192, 0.4687994511351 },
			{ 15.4211162384085, 12.8968882660397, 4.3267779995134, -1.8025500271447, 0.8151099343527 },
			{ 20.0518426536508, 13.4780491822807, 8.0545703111639, -1.4807768397938, 0.8330718096630 },
			{ 25.2756377203265, 14.0592100985217, 11.4058634577602, -0.1894358359554, 0.9988983039332 },
			{ 22.3399116516356, 14.5613210928137, 8.0505180814945, -0.2719275226726, 0.9949500881434 },
			{ 19.4330511346515, 15.0634320871057, 4.4146864348030, -0.0450673872573, 0.9986075895050 },
			{ 15.5743440275408, 15.5655430813978, -0.8309457281297, 0.8397466742728, 0.9513331929281 },
			{ 14.2152137310471, 16.1143236202088, -2.5032577528913, 0.6041478637295, 0.9802790804198 },
			{ 10.4466268193217, 16.6631041590199, -9.1111164128691, 2.8946390731709, 0.6361116722513 },
			{ 6.8363727508706, 17.2118846978310, -10.3253625568245, -0.0501493901358, 0.9998134288526 },
			{ 9.2248871573319, 17.7903433462106, -7.8626435829429, -0.7028126059358, 0.9918489900166 },
			{ 11.5170080038755, 18.3688019945901, -5.5711883348306, -1.2806056558840, 0.9036143700355 },
			{ 17.6596323618988, 18.9472606429697, 0.1248516164274, -1.4124798974983, 0.8396423896687 },
			{ 21.9815561700264, 19.4705628559132, 4.0721553326288, -1.5611620185156, 0.8267787563738 },
			{ 28.9000612233196, 19.9938650688567, 8.0548547142751, 0.8513414401877, 0.9679786505821 },
			{ 32.2440468292909, 20.5171672818002, 10.9856073304818, 0.7412722170089, 0.9652986697060 },
			{ 32.0124137923314, 21.0430516642624, 8.4103140318012, 2.5590480962678, 0.6913435115038 },
			{ 24.3632777472385, 21.5689360467246, 4.7802100670225, -1.9858683665086, 0.7609795420773 },
			{ 21.7892806474850, 22.0948204291868, -1.2317600178105, 0.9262202361087, 0.9404667128564 },
			{ 21.5716042249031, 22.5890578122113, -2.8679047358776, 1.8504511485694, 0.8185457042882 },
			{ 14.8817736687852, 23.0832951952358, -7.4224752454196, -0.7790462810310, 0.9515889001780 },
			{ 12.1617347215181, 23.5775325782603, -11.0346882717562, -0.3811095849860, 0.9942645306203 },
			{ 17.4956017010347, 23.9969747823523, -8.2822764353366, 1.7809033540189, 0.7872566380160 },
			{ 17.0797756823878, 24.4164169864444, -5.9375845542598, -1.3990567497967, 0.8923770320437 },
			{ 25.4356674141035, 24.8358591905364, -0.3359969137244, 0.9358051372915, 0.9596562394610 },
			{ 30.3052706879692, 25.2470311221379, 4.7624298351573, 0.2958097306740, 0.9990424083788 },
			{ 32.7685135659586, 25.6582030537394, 8.2107450919670, -1.1004345797478, 0.9340845959824 },
			{ 35.8602965626504, 26.0693749853409, 10.9793648352534, -1.1884432579440, 0.8933065786729 },
			{ 34.8002949236929, 26.5636238727732, 8.6734848229322, -0.4368137720125, 0.9794983866036 },
			{ 35.6396712902620, 27.0578727602055, 5.3020545671766, 3.2797439628799, 0.5625576684369 },
			{ 25.1017636136783, 27.5521216476379, -1.6272169107592, -0.8231411232004, 0.9611837345895 },
			{ 22.8180188403000, 28.1180572369798, -4.3187982869963, -0.9812401096835, 0.9337418671009 },
			{ 22.7315180720967, 28.6839928263218, -6.9524480774168, 0.9999733231918, 0.9565206405732 },
			{ 16.0932850984985, 29.2499284156638, -11.1198256253356, -2.0368176918297, 0.7768664914922 },
			{ 19.5586907877729, 29.8287080623234, -7.8096654488235, -2.4603518257270, 0.6930049701384 },
			{ 25.6416299983757, 30.4074877089830, -5.6701466788273, 0.9042889682200, 0.9607682159911 },
			{ 30.2628470809808, 30.9862673556426, -0.8261826088644, 0.1027623342027, 0.9997262240132 },
			{ 41.0388743777540, 31.6251683617577, 5.2692331643587, 4.1444728516376, 0.4199524331508 },
			{ 41.6343042052354, 32.2640693678728, 8.5396630395032, 0.8305717978593, 0.9642128857304 },
			{ 43.6799220667725, 32.9029703739880, 11.0078866016454, -0.2309349088608, 0.9836555441198 },
			{ 41.0850097483907, 33.5227579478548, 8.4554558177615, -0.8932040172256, 0.9443823488963 },
			{ 39.6750823915198, 34.1425455217217, 5.8159552848880, -0.2834184150899, 0.9824270599109 },
			{ 31.6122874004725, 34.7623330955886, -1.6498850297186, -1.5001606653975, 0.8743892420512 },
			{ 29.1621622497784, 35.3264886401658, -5.6022174011205, -0.5621089892669, 0.9742679975631 },
			{ 29.8318749562799, 35.8906441847430, -7.5815752233950, 1.5228059949319, 0.8862057540577 },
			{ 27.3169593005418, 36.4547997293201, -10.6399962772101, 1.5021558484318, 0.8735578702216 },
			{ 30.3803098165402, 37.0376895693819, -7.8266614136287, 1.1692816607870, 0.9197346315814 },
			{ 32.6898906858946, 37.6205794094438, -5.1049701150411, 0.1742813914920, 0.9995273107918 },
			{ 35.8775691758293, 38.2034692495056, -0.9795280369738, -1.3463720367025, 0.8889686354814 },
			{ 42.7972447632672, 38.7449409824910, 5.4144574828797, -1.3621537021035, 0.8109192911975 },
			{ 47.2380330372943, 39.2864127154765, 9.3062493021255, -1.3546289803077, 0.8845636110149 },
			{ 54.5214548424807, 39.8278844484620, 10.8609861754638, 3.8325842185549, 0.4033551204736 },
			{ 49.6705557480712, 40.3593441824895, 8.3270922520223, 0.9841193135594, 0.9565679050954 },
			{ 46.5368805320702, 40.8908039165171, 5.5165084314489, 0.1295681841042, 0.9999965297275 },
			{ 40.8793007016155, 41.4222636505447, -1.0839008592435, 0.5409379103142, 0.9875252295724 },
			{ 34.9747442848728, 41.9802441287924, -6.1721506675934, -0.8333491763262, 0.9528423578490 },
			{ 33.4863356980275, 42.5382246070402, -8.7541212967169, -0.2977676122958, 0.9943164414371 },
			{ 31.9907827865744, 43.0962050852880, -9.3208405425349, -1.7845817561787, 0.8391987269336 },
			{ 36.2541102958944, 43.6238270153547, -7.8867840880337, 0.5170673685734, 0.9758047869478 },
			{ 39.4985028703065, 44.1514489454214, -5.0402561369169, 0.3873100618020, 0.9870533026907 },
			{ 43.9819582396261, 44.6790708754881, -1.1721130766979, 0.4750004408358, 0.9745706290986 },
			{ 50.8231658435360, 45.2161020153221, 4.8999150983435, 0.7071487298704, 0.9732508938159 },
			{ 56.6876049996977, 45.7531331551561, 9.6288790812767, 1.3055927632649, 0.8742857332393 },
			{ 53.7044888571944, 46.2901642949900, 11.1434324280892, -3.7291078658848, 0.4205866068768 },
			{ 54.4121973474492, 46.8113091386879, 8.5790423556901, -0.9781541469287, 0.9598238692496 },
			{ 52.7460099929406, 47.3324539823857, 5.1386853348496, 0.2748706757052, 0.9847898564118 },
			{ 48.1369571932524, 47.8535988260836, -0.4669085875007, 0.7502669546695, 0.9558892256913 },
			{ 43.2760722384158, 48.3362291652832, -6.3171352726521, 1.2569783457848, 0.9007463042543 },
			{ 37.4371803202770, 48.8188595044827, -9.7845119031675, -1.5971672810382, 0.8788821518554 },
			{ 42.0270074809715, 49.3014898436823, -8.0208761163519, 0.7463937536411, 0.9574457319503 },
			{ 41.3627243185270, 49.7517447935070, -8.2414007401979, -0.1476197347821, 0.9994885863065 },
			{ 44.8896123684894, 50.2019997433317, -5.3230513483983, 0.0106639735561, 0.9999130488253 },
			{ 50.4640719308277, 50.6522546931564, -1.2374916996241, 1.0493089372955, 0.9262957667877 },
			{ 54.2211885693851, 51.0489447359274, 4.9378623146183, -1.7656184811606, 0.8331476489487 },
			{ 62.5525738041031, 51.4456347786985, 9.4491811545519, 1.6577578708527, 0.8279297206718 },
			{ 63.4253418894145, 51.8423248214696, 10.9334548687828, 0.6495621991621, 0.9558217129719 },
			{ 61.8577496655448, 52.2232409048262, 8.8261354100187, 0.8083733506999, 0.9651247801294 },
			{ 55.4889258192623, 52.6041569881827, 5.4242270383348, -2.5394582072552, 0.7015924617142 },
			{ 53.4830879341436, 52.9850730715393, -0.6719106489623, 1.1699255115667, 0.9232653948214 },
			{ 46.3439942426070, 53.3842025192624, -6.4072459905660, -0.6329622860893, 0.9823348261357 },
			{ 43.1553940257147, 53.7833319669855, -9.9410918475563, -0.6868460937145, 0.9779890700727 },
			{ 47.8625051601890, 54.1824614147085, -7.2522164062747, 0.9322601517552, 0.9401881162925 },
			{ 45.2605735419700, 54.6557583984503, -8.3704925843535, -1.0246922721268, 0.9369177925975 },
			{ 49.3869314069288, 55.1290553821921, -5.7413342951345, -0.0007896801288, 0.9996921898589 },
			{ 52.6173066266858, 55.6023523659338, -0.6746523133184, -2.3103934259296, 0.7386272221820 },
			{ 62.8569419952637, 56.1478023500594, 4.8392265392552, 1.8699131059492, 0.8094682281157 },
			{ 62.8957553956296, 56.6932523341849, 8.8296345774139, -2.6271315159691, 0.6751016860961 },
			{ 68.6688312235904, 57.2387023183104, 10.7378408147763, 0.6922880905037, 0.9515799657776 },
			{ 65.8616618077269, 57.7878694702290, 9.0994013304678, -1.0256089929700, 0.9338597641983 },
			{ 65.8795538039168, 58.3370366221476, 5.8342994548584, 1.7082177269108, 0.8258157608515 },
			{ 56.8141432151422, 58.8862037740662, -1.3402708463408, -0.7317897125832, 0.9677065134661 },
			{ 53.4689079768181, 59.4597171279688, -6.4295926363810, 0.4387834852303, 0.9837365569523 },
			{ 51.2976962551782, 60.0332304818714, -9.3732194212911, 0.6376851945979, 0.9814488822503 },
			{ 54.7021399639830, 60.6067438357740, -7.8172931432471, 1.9126892714561, 0.8211639609401 },
			{ 52.4769544857444, 61.1191499312501, -8.0340184612197, -0.6081769842859, 0.9529271901085 },
			{ 55.2603055496883, 61.6315560267261, -5.5429797933170, -0.8282706837208, 0.9497856193146 },
			{ 61.1230908925613, 62.1439621222022, -0.0393819614144, -0.9814892682265, 0.9560730971730 },
			{ 67.0950062779945, 62.6197349184303, 4.7832189603388, -0.3079476007746, 0.9958163929214 },
			{ 72.5235307815571, 63.0955077146585, 7.9922415549113, 1.4357815119873, 0.8607079217305 },
			{ 74.1658374011649, 63.5712805108866, 10.3928420836902, 0.2017148065881, 0.9951425607403 },
			{ 75.3946302291045, 64.1033675501705, 9.2008557064907, 2.0904069724432, 0.7836882077483 },
			{ 70.2944765528212, 64.6354545894543, 6.3625871205496, -0.7035651571828, 0.9719126941819 },
			{ 62.2989185185897, 65.1675416287382, -1.9985391602924, -0.8700839498561, 0.9521254803609 },
			{ 56.7262320247411, 65.7140918053384, -6.2658024055185, -2.7220573750788, 0.6487548102984 },
			{ 57.7351407889509, 66.2606419819385, -8.7813072436205, 0.2558060506329, 0.9993454382285 },
			{ 56.8212761847008, 66.8071921585387, -8.3787339465466, -1.6071820272913, 0.8246986751471 },
			{ 62.6604299795570, 67.3811754816226, -7.8559246917539, 3.1351791896883, 0.6394935052493 },
			{ 61.6451935638573, 67.9551588047065, -5.2176279431524, -1.0923372976969, 0.8977458121560 },
			{ 70.8589268303570, 68.5291421277905, 0.6411828513387, 1.6886018512279, 0.8074220082865 },
			{ 73.6787207050528, 69.1652525532692, 4.6932264485922, -0.1797582968086, 0.9999986938446 },
			{ 76.5623818680180, 69.8013629787480, 7.1424352510148, -0.3814163617448, 0.9994798313880 },
			{ 79.9054016284510, 70.4374734042268, 9.9596190282007, -0.4916908039764, 0.9880162985847 },
			{ 79.3873499272061, 71.0717659954274, 9.3134177860153, -0.9978338542366, 0.9475348469252 },
			{ 78.7410687642460, 71.7060585866281, 6.9097486861278, 0.1252614914901, 0.9979290381448 },
			{ 70.3301546008351, 72.3403511778287, -2.6788802398007, 0.6686836628071, 0.9668313065072 },
			{ 68.1785086337389, 72.9763599391604, -6.1025117637132, 1.3046604582917, 0.8567608388749 },
			{ 64.9662035129925, 73.6123687004921, -8.1770227106402, -0.4691424768594, 0.9868688732507 },
			{ 65.0445849452249, 74.2483774618237, -9.1084355898637, -0.0953569267351, 0.9991806127263 },
			{ 65.7903118858262, 74.8842576457032, -7.5873830049575, -1.5065627549194, 0.8091664150188 },
			{ 71.6003290025845, 75.5201378295826, -4.8318733281443, 0.9120645011462, 0.9628298432660 } };

}

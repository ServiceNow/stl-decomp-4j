package com.snc.ds.stats.stl;

import static java.lang.Math.PI;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CyclicSubSeriesSmootherTest {

	// Smoothing the cyclic sub-series extends the data one period in each direction. Ensure that when the data is
	// linear, that the extrapolations are linear.

	@Test
	public void TrendingSinusoidExtrapolationTest() {
		final int period = 24;
		double[] data = new double[2 * period];
		final double dx = 2 * PI / period;
		for (int i = 0; i < data.length; ++i) {
			final int amplitude = 10 - i / period;
			data[i] = amplitude * Math.sin(i * dx);
		}

		double[] extendedData = new double[4 * period];

		LoessSettings seasonalSettings = new LoessSettings(7); // Sub-cycle data is linear so width shouldn't matter

		CyclicSubSeriesSmoother sssmoother = new CyclicSubSeriesSmoother(seasonalSettings, data.length, period);

		sssmoother.smoothSeasonal(data, extendedData, null);

		for (int i = 0; i < extendedData.length; ++i) {
			final int amplitude = 11 - i / period; // An extra for the extrapolation before.
			final double value = amplitude * Math.sin(i * dx);
			assertEquals(String.format("test point %d", i), value, extendedData[i], 1.0e-11);
		}
	}

	@Test
	public void shouldExtrapolateFourPeriodsForwards() {
		final int period = 24;
		double[] data = new double[2 * period];
		final double dx = 2 * PI / period;
		for (int i = 0; i < data.length; ++i) {
			final int amplitude = 10 - i / period;
			data[i] = amplitude * Math.sin(i * dx);
		}

		double[] extendedData = new double[6 * period];

		LoessSettings seasonalSettings = new LoessSettings(7); // Sub-cycle data is linear so width shouldn't matter

		CyclicSubSeriesSmoother sssmoother = new CyclicSubSeriesSmoother(seasonalSettings, data.length, period, 4);
		sssmoother.smoothSeasonal(data, extendedData, null);

		for (int i = 0; i < extendedData.length; ++i) {
			final int amplitude = 10 - i / period;
			final double value = amplitude * Math.sin(i * dx);
			assertEquals(String.format("test point %d", i), value, extendedData[i], 1.0e-11);
		}
	}

	@Test
	public void shouldExtrapolateTwoPeriodsBackwardAndTwoPeriodsForward() {
		final int period = 24;
		double[] data = new double[2 * period];
		final double dx = 2 * PI / period;
		for (int i = 0; i < data.length; ++i) {
			final int amplitude = 10 - i / period;
			data[i] = amplitude * Math.sin(i * dx);
		}

		double[] extendedData = new double[6 * period];

		LoessSettings seasonalSettings = new LoessSettings(7); // Sub-cycle data is linear so width shouldn't matter

		CyclicSubSeriesSmoother sssmoother = new CyclicSubSeriesSmoother(seasonalSettings, data.length, period, 2, 2);
		sssmoother.smoothSeasonal(data, extendedData, null);

		for (int i = 0; i < extendedData.length; ++i) {
			final int amplitude = 12 - i / period; // Two extra for the extrapolation before.
			final double value = amplitude * Math.sin(i * dx);
			assertEquals(String.format("test point %d", i), value, extendedData[i], 1.0e-11);
		}
	}

}

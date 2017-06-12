package com.github.servicenow.ds.stats.stl;

import org.junit.Test;

import static com.github.servicenow.ds.stats.stl.CyclicSubSeriesSmoother.Builder;
import static java.lang.Math.PI;
import static org.junit.Assert.assertEquals;

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

		Builder builder = new Builder();
		builder.setWidth(7); // Sub-cycle data is linear so width shouldn't matter

		CyclicSubSeriesSmoother sssmoother = builder.setDataLength(data.length).setPeriodicity(period)
				.extrapolateForwardAndBack(1).build();

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

		Builder builder = new Builder();
		builder = builder.setWidth(7); // Sub-cycle data is linear so width shouldn't matter
		builder = builder.extrapolateForwardOnly(4);

		CyclicSubSeriesSmoother sssmoother = builder.setDataLength(data.length).setPeriodicity(period).build();

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

		Builder builder = new Builder();
		builder = builder.setWidth(7); // Sub-cycle data is linear so width shouldn't matter
		builder = builder.setNumPeriodsForward(2).setNumPeriodsBackward(2);

		CyclicSubSeriesSmoother sssmoother = builder.setDataLength(data.length).setPeriodicity(period).build();

		sssmoother.smoothSeasonal(data, extendedData, null);

		for (int i = 0; i < extendedData.length; ++i) {
			final int amplitude = 12 - i / period; // Two extra for the extrapolation before.
			final double value = amplitude * Math.sin(i * dx);
			assertEquals(String.format("test point %d", i), value, extendedData[i], 1.0e-11);
		}
	}

	@Test(expected = IllegalArgumentException.class)
	public void degreeMustBePositive() {
		Builder builder = new Builder();
		builder.setDegree(-1);
	}

	@Test(expected = IllegalArgumentException.class)
	public void degreeMustBeLessThanThree() {
		Builder builder = new Builder();
		builder.setDegree(3);
	}

	@Test(expected = IllegalArgumentException.class)
	public void widthMustBeSet() {
		Builder builder = new Builder();
		builder.setDataLength(100).extrapolateForwardAndBack(1).setPeriodicity(12).build();
	}

	@Test(expected = IllegalArgumentException.class)
	public void dataLengthMustBeSet() {
		Builder builder = new Builder();
		builder.setWidth(3).extrapolateForwardAndBack(1).setPeriodicity(12).build();
	}

	@Test(expected = IllegalArgumentException.class)
	public void periodMustBeSet() {
		Builder builder = new Builder();
		builder.setDataLength(100).extrapolateForwardAndBack(1).setWidth(11).build();
	}

	@Test(expected = IllegalArgumentException.class)
	public void backwardExtrapolationMustBeSet() {
		Builder builder = new Builder();
		builder.setDataLength(100).setNumPeriodsForward(1).setWidth(11).setPeriodicity(12).build();
	}

	@Test(expected = IllegalArgumentException.class)
	public void forwardExtrapolationMustBeSet() {
		Builder builder = new Builder();
		builder.setDataLength(100).setNumPeriodsBackward(1).setWidth(11).setPeriodicity(12).build();
	}
}

package com.github.servicenow.ds.stats;

import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Tests for TimeSeriesUtilities
 *
 * Created by Jim Crotinger on 20-Apr-2016.
 */

public class TimeSeriesUtilitiesTest {
    @Test
    public void smaWithWindowEqualLengthIsJustAverage() {
        int length = (int) (Math.random() * 1000 + 1); // uniform random in [1..1000]

        double[] data = createRandomArray(length);

        double sum = 0.0;
        for (int i = 0; i < data.length; ++i) {
            sum += data[i];
        }

        double mean = sum / data.length;

        double[] average = TimeSeriesUtilities.simpleMovingAverage(data, data.length);

        assertEquals("average has length of 1", 1, average.length);
        assertEquals("average[0] value is just the average", mean, average[0], 1.0e-10);
    }

    @Test
    public void smaWithWindowEqualOneIsJustData() {
        double[] data = createRandomArray(10);

        double[] average = TimeSeriesUtilities.simpleMovingAverage(data, 1);

        assertEquals("average has length of data.length", data.length, average.length);
        for (int i = 0; i < data.length; ++i) {
            assertEquals("average is just the original data", data[i], average[i], 1.0e-10);
        }
    }

    @Test
    public void smaRandomDataTest() {
        int length = (int) (Math.random() * 1000 + 1); // uniform random in [1..1000]

        double[] data = createRandomArray(length);

        int window = (int) (Math.random() * length);
        window = Math.max(window, 2);
        window = Math.min(window, length);

        double[] average = TimeSeriesUtilities.simpleMovingAverage(data, window);

        assertEquals("average has right length", data.length - window + 1, average.length);

        for (int i = 0; i < average.length; ++i) {
            double sum = 0.0;
            for (int j = 0; j < window; ++j) {
                sum += data[i + j];
            }
            double mean = sum / window;

            assertEquals("moving average is correct", mean, average[i], 1.0e-10);
        }
    }

    @Test(expected=IllegalArgumentException.class)
    public void lengthConsistencyTest() {
        int length = (int) (Math.random() * 1000 + 1); // uniform random in [1..1000]

        double[] data = createRandomArray(length);

        int window = (int) (Math.random() * length);
        window = Math.max(window, 2);
        window = Math.min(window, length);

        double[] average = new double[data.length - window - 1]; // Sign mistake in size...
        TimeSeriesUtilities.simpleMovingAverage(data, window, average);
    }

    private double[] createRandomArray(int length) {
        double[] data = new double[length];
        for (int i = 0; i < data.length; ++i) {
            data[i] = Math.random() * 100 - 50; // [-50..50];
        }
        return data;
    }

}
package com.github.servicenow.ds.stats;

/**
 * TimeSeriesUtilities is a collection of static functions implementing common time-series utilities.
 *
 * Created by Jim Crotinger on 19-Apr-2016.
 */
public class TimeSeriesUtilities {

    /**
     * Compute the simple moving average of the array data using the specified window size and return the results.
     *
     * @param data double[] array of input data
     * @param window int width of moving average
     * @return double[] of length data.length - window + 1, containing the simple moving average
     */
    public static double[] simpleMovingAverage(final double[] data, final int window) {
        double[] average = new double[data.length - window + 1];
        simpleMovingAverage(data, window, average);
        return average;
    }

    /**
     * Compute the simple moving average of the array data using the specified window size and return the results
     * in the array average.
     *
     * @param data double[] array of input data
     * @param window int width of moving average
     * @param average double[] array of output data of length at least data.length - window + 1
     */
    public static void simpleMovingAverage(final double[] data, final int window, final double[] average) {
		if (average.length < data.length - window + 1)
		    throw new IllegalArgumentException("simpleMovingAverage: insufficient memory to store moving average");

        // The simple moving average picks up one point from the first "window" points of the original data and then
        // data.length - window additional points for the rest of the data.

        double windowSum = 0.0;
        for (int i = 0; i < window; ++i)
            windowSum += data[i];

        average[0] = windowSum / window;

        // Now roll through the additional data subtracting the contribution from the index that has left the window
        // and adding the contribution from the next index to enter the window. Last window above was [0, window - 1],
        // so we need to start by removing data[0] and adding data[window] and move forward until we add the
        // last data point; i.e. until windowEnd == data.length - 1.
        int windowEnd = window;
        int windowStart = 0;
        for (int j = 1; j < data.length - window + 1; ++j) {
            // loops data.length - window + 1 - 1 = data.length - window times
            windowSum += (data[windowEnd] - data[windowStart]);
            ++windowStart;
            ++windowEnd;
            average[j] = windowSum / window;
        }
    }
}

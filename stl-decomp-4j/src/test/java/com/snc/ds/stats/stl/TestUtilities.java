package com.snc.ds.stats.stl;

import java.io.FileWriter;
import java.io.IOException;

/**
 * TestUtilities - exactly that.
 *
 * Created by Jim Crotinger on 2-Feb-17.
 */
class TestUtilities {

	static void dumpStlResultsToFile(
			double[] data, SeasonalTrendLoess.Decomposition stl, String fileName) throws IOException {
		double[] trend = stl.getTrend();
		double[] seasonal = stl.getSeasonal();
		double[] residuals = stl.getResiduals();

		FileWriter out = new FileWriter(fileName);
		for (int i = 0; i < data.length; ++i) {
			out.write(String.format("%17.14f, %17.14f, %17.14f, %17.14f\n", data[i], seasonal[i], trend[i],
					residuals[i]));
		}
		out.close();
	}

}

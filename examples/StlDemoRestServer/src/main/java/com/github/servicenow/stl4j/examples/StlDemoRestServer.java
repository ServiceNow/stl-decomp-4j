package com.github.servicenow.stl4j.examples;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import com.fasterxml.jackson.databind.ObjectMapper;

import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;

import com.github.servicenow.ds.stats.stl.SeasonalTrendLoess;

import static spark.Spark.*;

/**
 * StlDemoRestServer - simple Spark REST server to serve up STL results to a web API
 *
 * Created by Jim Crotinger on 24-Mar-2017
 */
public class StlDemoRestServer {

	public static void main(String[] args) throws IOException {

		StlResults results = new StlResults();

		TimeSeries ts = getTimeSeries("co2.csv");

		SeasonalTrendLoess.Builder stlBuilder = new SeasonalTrendLoess.Builder();

		final double[] vs = ts.values.stream().mapToDouble(i -> i).toArray();

		for (int i = 0; i < ts.values.size(); ++i)
			vs[i] = ts.values.get(i);

		final SeasonalTrendLoess smoother = stlBuilder.
				setPeriodLength(12).
				setSeasonalWidth(35).
				setNonRobust().
				buildSmoother(vs);

		final SeasonalTrendLoess.Decomposition stl = smoother.decompose();

//		printResults(ts, stl);

		results.time = ts.times.stream().mapToLong(i -> i).toArray();
		results.value = vs;
		results.seasonal = stl.getSeasonal();
		results.trend = stl.getTrend();
		results.residual = stl.getResidual();
		results.weight = stl.getWeights();

		ObjectMapper mapper = new ObjectMapper();

		String results_json = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(results);

		enableCORS("*", "*", "*");

		get("/stldemo", (req, res) -> results_json);

		System.out.println("Done");
	}

	private static void printResults(TimeSeries ts, SeasonalTrendLoess.Decomposition stl) {
		for (int i = 0; i < ts.values.size(); ++i) {
			System.out.println(
					String.format("time = %d, value = %f, seasonal = %f, trend = %f, residual = %f, weights = %f",
							ts.times.get(i), ts.values.get(i),
							stl.getSeasonal()[i], stl.getTrend()[i], stl.getResidual()[i], stl.getWeights()[i]));
		}
	}

	public static class TimeSeries {
		ArrayList<Double> values = new ArrayList<>();
		ArrayList<Long> times = new ArrayList<>();
	}

	@SuppressWarnings("Duplicates")
	public static TimeSeries getTimeSeries(String fileName) throws IOException {
		CSVReaderBuilder builder = new CSVReaderBuilder(new FileReader(fileName));

		TimeSeries ts = new TimeSeries();

		try (CSVReader reader = builder.withSkipLines(1).build()) {

			String[] nextLine;
			while ((nextLine = reader.readNext()) != null) {
				double dateAsYear = Double.parseDouble(nextLine[1]);
				long time = (long) ((dateAsYear - 1970.0) * 365.25 * 24 * 60 * 60 * 1000);
				ts.times.add(time);

				double value = Double.parseDouble(nextLine[2]);
				ts.values.add(value);
			}
		}
		return ts;
	}

	// Enables CORS on requests. This method is an initialization method and should be called once.
	private static void enableCORS(final String origin, final String methods, final String headers) {

		options("/*", (request, response) -> {

			String accessControlRequestHeaders = request.headers("Access-Control-Request-Headers");
			if (accessControlRequestHeaders != null) {
				response.header("Access-Control-Allow-Headers", accessControlRequestHeaders);
			}

			String accessControlRequestMethod = request.headers("Access-Control-Request-Method");
			if (accessControlRequestMethod != null) {
				response.header("Access-Control-Allow-Methods", accessControlRequestMethod);
			}

			return "OK";
		});

		before((request, response) -> {
			response.header("Access-Control-Allow-Origin", origin);
			response.header("Access-Control-Request-Method", methods);
			response.header("Access-Control-Allow-Headers", headers);
			// Note: this may or may not be necessary in your particular application
			response.type("application/json");
		});
	}
}
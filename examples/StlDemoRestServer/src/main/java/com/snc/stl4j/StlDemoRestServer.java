package com.snc.stl4j;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.codehaus.jackson.map.ObjectMapper;

import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.snc.ds.stats.stl.SeasonalTrendLoess;

import static spark.Spark.*;

/**
 * StlDemoRestServer - simple Spark REST server to serve up STL results to a web API
 *
 * Created by Jim Crotinger on 24-Mar-2017
 */
public class StlDemoRestServer {
	public static void main(String[] args) throws IOException {

		StlResults results = new StlResults();

		CSVReaderBuilder builder = new CSVReaderBuilder(new FileReader("co2.csv"));

		try (CSVReader reader = builder.withSkipLines(1).build()) {

			ArrayList<Double> values = new ArrayList<>();
			ArrayList<Long> times = new ArrayList<>();

			String [] nextLine;
			while ((nextLine = reader.readNext()) != null) {
				double dateAsYear = Double.parseDouble(nextLine[1]);
				long time = (long) ((dateAsYear - 1970.0) * 365.25 * 24 * 60 * 60 * 1000);
				times.add(time);

				double co2 = Double.parseDouble(nextLine[2]);
				values.add(co2);
			}

			SeasonalTrendLoess.Builder stlBuilder = new SeasonalTrendLoess.Builder();

			final double[] vs = values.stream().mapToDouble(i->i).toArray();

			for (int i = 0; i < values.size(); ++i)
				vs[i] = values.get(i);

			final SeasonalTrendLoess smoother = stlBuilder.
					setPeriodLength(12).
					setSeasonalWidth(35).
					setNonRobust().
					buildSmoother(vs);

			final SeasonalTrendLoess.Decomposition stl = smoother.decompose();

//			for (int i = 0; i < values.size(); ++i) {
//				System.out.println(
//						String.format("time = %d, value = %f, seasonal = %f, trend = %f, residual = %f, weights = %f",
//								times.get(i), values.get(i),
//								stl.getSeasonal()[i], stl.getTrend()[i], stl.getResiduals()[i], stl.getWeights()[i]));
//			}

			results.time = times.stream().mapToLong(i->i).toArray();
			results.value = vs;
			results.seasonal = stl.getSeasonal();
			results.trend = stl.getTrend();
			results.residual = stl.getResiduals();
			results.weight = stl.getWeights();
		}

		ObjectMapper mapper = new ObjectMapper();

		String results_json = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(results);

		enableCORS("*", "*", "*");

		get("/stldemo", (req, res) -> results_json);

		System.out.println("Done");
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
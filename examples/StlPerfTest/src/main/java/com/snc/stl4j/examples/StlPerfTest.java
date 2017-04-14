package com.snc.stl4j.examples;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.DefaultParser;

import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.snc.ds.stats.stl.SeasonalTrendLoess;

public class StlPerfTest {
	private static int fTimedIterations = 30;
	private static int fWarmupIterations = 2000;
	private static String fDataFilePath = "../StlDemoRestServer/co2.csv";

	public static void main(String[] args) throws IOException {
		parseCommandLine(args);

		TimeSeries ts = getTimeSeries(fDataFilePath);

		SeasonalTrendLoess.Builder stlBuilder = new SeasonalTrendLoess.Builder();

		final double[] vs = ts.values.stream().mapToDouble(i -> i).toArray();

		for (int i = 0; i < ts.values.size(); ++i)
			vs[i] = ts.values.get(i);

		final SeasonalTrendLoess smoother = stlBuilder.
				setPeriodLength(12).
				setSeasonalWidth(35).
				setNonRobust().
				buildSmoother(vs);

		for (int i = 0; i < fWarmupIterations; ++i) {
			//noinspection unused
			SeasonalTrendLoess.Decomposition unused = smoother.decompose();
		}

		long start = System.nanoTime();
		for (int i = 0; i < fTimedIterations; ++i) {
			//noinspection unused
			SeasonalTrendLoess.Decomposition unused = smoother.decompose();
		}

		double elapsed = (System.nanoTime() - start) * 1.0e-9;
		System.out.println(
				String.format("Elapsed time = %f s; Time per iteration = %f ms",
						elapsed, 1000 * elapsed / fTimedIterations));
	}

	private static class TimeSeries {
		final ArrayList<Double> values = new ArrayList<>();
		final ArrayList<Long> times = new ArrayList<>();
	}

	@SuppressWarnings("Duplicates")
	private static TimeSeries getTimeSeries(String fileName) throws IOException {
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


	private static void parseCommandLine(String[] args) {
		Options options = new Options();

		Option input = new Option("n", "timed-iterations", true, "number of iterations of timing loop");
		input.setRequired(false);
		options.addOption(input);

		Option output = new Option("w", "warmup-iterations", true, "number of warm-up iterations before timing loop");
		output.setRequired(false);
		options.addOption(output);

		Option filename = new Option("f", "filename", true, "data file to use for timing");
		filename.setRequired(false);
		options.addOption(filename);

		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();
		CommandLine cmd;

		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			formatter.printHelp("StlPerfTest", options);

			System.exit(1);
			return;
		}

		String nStr = cmd.getOptionValue("number");
		if (nStr != null)
			fTimedIterations = Integer.parseInt(nStr);

		String wStr = cmd.getOptionValue("warmup-iterations");
		if (wStr != null)
			fWarmupIterations = Integer.parseInt(wStr);

		String fStr = cmd.getOptionValue("filename");
		if (fStr != null)
			fDataFilePath = fStr;
	}
}

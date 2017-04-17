package com.snc.stl4j.examples;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.snc.ds.stats.stl.SeasonalTrendLoess;

public class StlPerfTest {
	private static int fTimedIterations;
	private static int fWarmupIterations;
	private static boolean fRunCo2 = true;

	public static void main(String[] args) throws IOException {
		parseCommandLine(args);

		if (fRunCo2)
			runCo2Test();
		else
			runHourlyTest();
	}

	private static void runCo2Test() throws IOException {
		TimeSeries ts = getCo2Data();

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

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		dumpStlDecomposition(stl);
	}

	private static void runHourlyTest() throws IOException {
		TimeSeries ts = getHourlyData();

		SeasonalTrendLoess.Builder stlBuilder = new SeasonalTrendLoess.Builder();

		final double[] vs = ts.values.stream().mapToDouble(i -> i).toArray();

		for (int i = 0; i < ts.values.size(); ++i)
			vs[i] = ts.values.get(i);

		final SeasonalTrendLoess smoother = stlBuilder.
				setPeriodLength(8736).
				setSeasonalWidth(893451).
				setSeasonalDegree(0).
				setNonRobust().
				setTrendWidth(13105). // defaults might work for these
				setLowpassWidth(8737).
				setSeasonalJump(89346).
				setTrendJump(1311).
				setLowpassJump(874).
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

		SeasonalTrendLoess.Decomposition stl = smoother.decompose();

		dumpStlDecomposition(stl);
	}

	private static void dumpStlDecomposition(SeasonalTrendLoess.Decomposition stl)
			throws FileNotFoundException, UnsupportedEncodingException {
		try (PrintWriter writer = new PrintWriter("output.csv", "UTF-8")) {
			for (int i = 0; i < stl.getData().length; ++i) {
				double d = stl.getData()[i];
				double s = stl.getSeasonal()[i];
				double t = stl.getTrend()[i];
				double r = stl.getResidual()[i];
				writer.println(String.format("%.17E, %.17E, %.17E, %.17E", d, s, t, r));
			}
		}
	}

	private static class TimeSeries {
		final ArrayList<Double> values = new ArrayList<>();
		final ArrayList<Long> times = new ArrayList<>();
	}

	@SuppressWarnings("Duplicates")
	private static TimeSeries getCo2Data() throws IOException {

		final String path =  "../StlDemoRestServer/co2.csv";

		CSVReaderBuilder builder = new CSVReaderBuilder(new FileReader(path));

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

	private static TimeSeries getHourlyData() throws IOException {
		final String path = "./fortran_benchmark/hourly_stl_test.csv";

		CSVReaderBuilder builder = new CSVReaderBuilder(new FileReader(path));

		TimeSeries ts = new TimeSeries();

		try (CSVReader reader = builder.build()) {

			String[] nextLine;
			long time = 1492457959000L;
			while ((nextLine = reader.readNext()) != null) {
				ts.times.add(time);
				time += 3600 * 1000;
				double value = Double.parseDouble(nextLine[0]);
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

		Option hourly = new Option("h", "hourly", false, "whether to use hourly data");
		hourly.setRequired(false);
		options.addOption(hourly);

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

		if (cmd.hasOption("hourly")) {
			System.out.println("Running hourly stress test");
			fRunCo2 = false;
			fTimedIterations = 200;
			fWarmupIterations = 30;
		} else {
			System.out.println("Running CO2 test");
			fTimedIterations = 2000;
			fWarmupIterations = 30;
		}

		String nStr = cmd.getOptionValue("number");
		if (nStr != null)
			fTimedIterations = Integer.parseInt(nStr);

		String wStr = cmd.getOptionValue("warmup-iterations");
		if (wStr != null)
			fWarmupIterations = Integer.parseInt(wStr);

	}
}

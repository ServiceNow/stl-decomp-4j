# Seasonal Decomposition of Time Series by Loess

The Seasonal-Trend-Loess (STL) algorithm decomposes a time series into seasonal, trend and residual components. The algorithm uses [Loess interpolation](https://en.wikipedia.org/wiki/Local_regression) (original paper [here](http://www.stat.washington.edu/courses/stat527/s13/readings/Cleveland_JASA_1979.pdf)) to smooth the cyclic sub-series (e.g. all January values in the CO<sub>2</sub> data). After removing the seasonality from the signal, the remainder is smoothed (in multiple steps) to find the trend. This process is repeated and may include robustness iterations that take advantage of the weighted-least-squares underpinnings of Loess to remove the effects of outliers. The details are described in [STL: A Seasonal-Trend Decomposition Procedure Based on Loess](http://www.wessa.net/download/stl.pdf).   

**_stl-decomp-4j_** is a Java port of the original Ratfor/Fortran available from [Netlib](http://netlib.org) ([original source here](http://netlib.org/a/stl)), extended to support local quadratic interpolation.

As with the original Fortran version (and the [R](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/stl.html) and [Python](https://github.com/jcrotinger/pyloess) versions, which use the Fortran version under the hood), **_stl-decomp-4j_** expects equally spaced data with no missing values.

## Example
```java
double[] values = getSomeMonthlyData(); // Monthly time-series data

SeasonalTrendLoess.Builder builder = new SeasonalTrendLoess.Builder();
SeasonalTrendLoess smoother = builder.
			setPeriodLength(12).    // Data has a period of 12
			setSeasonalWidth(35).   // Monthly data smoothed over 35 years
			setNonRobust().         // Not expecting outliers, so no robustness iterations
			buildSmoother(values);

SeasonalTrendLoess.Decomposition stl = smoother.decompose();
double[] seasonal = stl.getSeasonal();
double[] trend = stl.getTrend();
double[] residual = stl.getResidual();
```

The `examples/StlDemoRestServer` directory includes a copy of the [Monthly CO<sub>2</sub> Measurement Data](http://www.esrl.noaa.gov/gmd/ccgg/trends/) and a simple REST server that reads this data, performs an STL decomposition on the data, and serves up the results to `http://localhost:4567/stldemo`. The `examples/StlDemoRestServer/index.html` file loads the data from this server and plots the resulting decomposition.

##Documentation
JavaDoc

## Interface

At a minimum, the STL algorithm requires specifying the periodicity of the data (e.g. 12 for monthly) and the width of the Loess smoother used to smooth the cyclic seasonal sub-series. In general, there are three Loess smoothers that each require three parameters --- a width, a degree, and a jump. The width specifies the number of data points that the local interpolation uses to smooth each point, the degree specifies the degree of the local polynomial that is fit to the data, and the jump specifies how many points are skipped between Loess interpolations, with linear interpolation being done between these points. Because of this complexity, construction is done via Builder objects, as shown in the simple example above.

![CO2 Plot](examples/StlDemoRestServer/co2_stl_highchart.jpg)

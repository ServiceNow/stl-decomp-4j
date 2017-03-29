# stl-decomp-4j

stl-decomp-4j is a Java implementation of the seasonal-trend decomposition algorith described in [STL: A Seasonal-Trend Decomposition Procedure Based on Loess](http://www.wessa.net/download/stl.pdf). This version is a Java port of the original Ratfor/Fortran available from [Netlib](http://netlib.org/a/stl), with an extension to support local quadratic interpolation [the underlying LOESS interpolation can do flat (degree = 0), linear (degree = 1) or quadratic (degree = 2) local interpolation].

As with the original Fortran version, this version of the STL algorithm expects equally spaced data with no missing values.

## Example

```
double[] values; // ...

SeasonalTrendLoess smoother = stlBuilder.
    setPeriodLength(12).
    setSeasonalWidth(35).
    setNonRobust().
    buildSmoother(values);

SeasonalTrendLoess.Decomposition stl = smoother.decompose();
```

The `examples/StlDemoRestServer` directory includes a copy of the [Monthly CO_2 Measurement Data](http://www.esrl.noaa.gov/gmd/ccgg/trends/) and a simple REST server that reads this data, performs an STL decomposition on the data, and serves up the results to http://localhost:4567/stldemo. The `index.html` file in that directory will load data from this endpoint and plot the resulting decomposition.

![CO2 Plot](examples/StlDemoRestServer/co2_stl_highchart.jpg)

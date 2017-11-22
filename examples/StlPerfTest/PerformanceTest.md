# Performance Test

## Java Tests

* Build the Java version with `mvn clean install`.
* Run the Java CO<sub>2</sub> test with `mvn -e test -Pco2`
* Run the Java Hourly test with `mvn -e test -Phourly`

## Fortran Tests

From the `fortran_benchmark` subdirectory, do

* `make run-benchmark` to build and run the CO<sub>2</sub> test
* `make run-hourly` to build and run the Hourly test

## Jupyter Notebook

The notebook `StlJavaFortranComparison.ipynb` compares the results from Java and Fortran STL implementation on the CO<sub>2</sub> data. It can be run by running

```
$ jupyter notebook . 
```
in this directory.

Note that the hourly benchmarks overwrite the test data with hourly data, so be sure to re-run the Java and Fortran CO<sub>2</sub> examples after running the hourly example (or modify the notebook accordingly).

# Performance Test

## Java Tests

* Build the Java version with `mvn clean install`.
* Run the Java CO<sub>2</sub> test with `mvn -e test -Pco2`
* Run the Java Hourly test with `mvn -e test -Phourly`

## Fortran Tests

From the `fortran_benchmark` subdirectory, do

* `make run-benchmark` to build and run the CO<sub>2</sub> test
* `make run-hourly` to build and run the Hourly test



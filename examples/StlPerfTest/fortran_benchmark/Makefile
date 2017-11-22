O= sample.o stl.o 
FC=gfortran
LIBDIR=/usr/lib
FFLAGS=-O3 -fdefault-real-8

clean:
	rm -f *.o a.out benchmark hourly

out: a.out co2
	./a.out < co2 >stl.txt

a.out: $O
	$(FC) -O3 -L/usr/lib $O

benchmark: benchmark.o stl.o
	$(FC) -O3 -L/usr/lib -o benchmark benchmark.o stl.o

run-benchmark: benchmark
	./benchmark <co2_2017

hourly: hourly.o stl.o
	$(FC) -O3 -L/usr/lib -o hourly hourly.o stl.o

run-hourly: hourly
	./hourly <hourly_stl_test.csv



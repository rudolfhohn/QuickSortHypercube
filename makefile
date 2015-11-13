# P = number of processors
P=4 

hyperquick: hyperquick.cc
	mpiCC $^ -o $@

run: hyperquick
	mpirun -np $(P) $^ 17 100

clean:
	rm -f *.o

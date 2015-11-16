# P = number of processors
P=2 

hyperquick: hyperquick.cc
	mpiCC $^ -o $@

run: hyperquick
	mpirun -np $(P) $^ 17 20

clean:
	rm -f *.o

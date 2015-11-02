# P = number of processors
P=4 

hyperquick.o: hyperquick.cc
	mpiCC $^ -o $@

run: hyperquick.o
	mpirun -np $(P) $^ 17 5

clean:
	rm -f *.o

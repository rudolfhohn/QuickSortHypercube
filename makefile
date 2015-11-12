# P = number of processors
P=4 

hyperquick.o: hyperquick.cc
	mpiCC $^ -o $@

run: hyperquick.o
	mpirun -np $(P) $^ 17 $(P)

clean:
	rm -f *.o

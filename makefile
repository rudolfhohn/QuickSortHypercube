# P = number of processors
P=4 

hyperquick.o: hyperquick.cc
	mpiCC $^ -o $@

run: hyperquick.o
	mpirun -np $(P) $^

clean:
	rm -f *.o

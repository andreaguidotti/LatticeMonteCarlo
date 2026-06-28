CC=gcc
CFLAGS=-O3 -Wall -Wextra   -Wconversion 

binaries= cluster random geometry obs clusterSim bootstrap
all: $(binaries)

random.o: ../lib/random.c  ../include/random.h
	$(CC) $(CFLAGS) -c ../lib/random.c

geometry.o: ../lib/geometry.c  ../include/geometry.h
	$(CC) $(CFLAGS) -c ../lib/geometry.c	

cluster: ../src/isingCluster.c random.o geometry.o
	$(CC) $(CFLAGS) -c ../src/isingCluster.c
	$(CC) $(CFLAGS) isingCluster.o random.o geometry.o -o $@ -lm



obs: ../src/observables.c
	$(CC) $(CFLAGS) ../src/observables.c -o $@ -lm 

.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o

CC=gcc
CFLAGS=-O3 -Wall -Wextra   -Wconversion 

binaries= z2model obs

all: $(binaries)

random.o: ../lib/random.c  ../include/random.h
	$(CC) $(CFLAGS) -c ../lib/random.c

geometry.o: ../lib/geometry.c  ../include/geometry.h
	$(CC) $(CFLAGS) -c ../lib/geometry.c	

z2model: ../src/z2model.c random.o geometry.o
	$(CC) $(CFLAGS) -c ../src/z2model.c
	$(CC) $(CFLAGS) z2model.o random.o geometry.o -o $@ -lm

obs: ../src/observables.c 
	$(CC) $(CFLAGS) -c ../src/observables.c
	$(CC) $(CFLAGS) observables.o -o $@ -lm	

.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o

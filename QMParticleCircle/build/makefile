CC=gcc
CFLAGS=-O3 -Wall -Wextra   -Wconversion 

binaries=  obs windingMC
all: $(binaries)

random.o: ../lib/random.c  ../include/random.h
	$(CC) $(CFLAGS) -c ../lib/random.c

geometry.o: ../lib/geometry.c  ../include/geometry.h
	$(CC) $(CFLAGS) -c ../lib/geometry.c

windingMC: ../src/windingMC.c random.o 
	$(CC) $(CFLAGS) -c ../src/windingMC.c
	$(CC) $(CFLAGS)windingMC.o random.o  -o $@ -lm	

obs: ../src/observable.c
	$(CC) $(CFLAGS) ../src/observable.c -o $@ -lm 

.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o

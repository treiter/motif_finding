CC=gcc

all: motif_finder.c
	gcc -Wall -o motif_finder motif_finder.c

clean:
	rm motif_finder
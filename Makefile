CC=gcc
CFLAGS= -std=c99 -Wall -lm

all: motif_finder.c
	$(CC) $(CFLAGS) -o motif_finder motif_finder.c

clean:
	rm motif_finder
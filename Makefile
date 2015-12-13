CC=gcc
CFLAGS= -std=c99 -Wall -lm -DTRACE -D_BSD_SOURCE -O3

all: motif_finder.c
	$(CC) $(CFLAGS) -o MotifFinder motif_finder.c

clean:
	rm MotifFinder
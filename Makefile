CC = gcc
FLAGS = -O2

.PHONY: all

all: qwest

qwest: qwest.o
	$(CC) -o $@ $<

%.o: %.c
	$(CC) -c -o $@ $< $(FLAGS)

clean: 
	rm -f qwest *.o

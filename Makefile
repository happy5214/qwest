CC = gcc
FLAGS = -O2 -Wall

objs = qwest.o int128.o carg_parser.o

.PHONY: all clean

all: qwest

qwest: $(objs)
	$(CC) -o $@ $^

%.o: %.c
	$(CC) -c -o $@ $< $(FLAGS)

clean: 
	rm -f qwest *.o

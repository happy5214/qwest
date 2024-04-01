CC = gcc
FLAGS = -O2 -Wall

objs = qwest.o int128.o

.PHONY: all

all: qwest

qwest: $(objs)
	$(CC) -o $@ $(objs)

%.o: %.c
	$(CC) -c -o $@ $< $(FLAGS)

clean: 
	rm -f qwest *.o

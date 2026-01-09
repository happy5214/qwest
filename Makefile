CPP = g++
FLAGS = -O2 -Wall

qwest_objs = qwest.o int128.o
arg_parser_objs = arg_parser.o

qwest_headers = int128.h
arg_parser_headers = arg_parser.h

.PHONY: all clean

all: qwest

qwest: $(qwest_objs) $(arg_parser_objs)
	$(CPP) -o $@ $^

%.o: %.cpp $(qwest_headers) $(arg_parser_headers)
	$(CPP) $(FLAGS) -c $< -o $@

%.o: %.cc $(arg_parser_headers)
	$(CPP) $(FLAGS) -c $< -o $@

clean: 
	rm -f qwest *.o

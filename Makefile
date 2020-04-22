all: (PROG) (PROG2)

PROG = main
CC = g++
FLAGS = -std=c++11 -Wall -g
OBJS = main.o viterbi.o parallel_viterbi.o util.o

PROG2 = generator

(PROG) : $(OBJS)
	$(CC) $(FLAGS) -o $(PROG) $(OBJS)
	rm -f $(OBJS)

(PROG2) : $(OBJS2)
	$(CC) $(FLAGS) -o $(PROG2) $(OBJS2)
	rm -f $(OBJS2)

main.o :
	$(CC) $(FLAGS) -c main.cpp

viterbi.o : viterbi.h
	$(CC) $(FLAGS) -c viterbi.cpp

parallel_viterbi.o : parallel_viterbi.h
	$(CC) $(FLAGS) -c parallel_viterbi.cpp

util.o : util.h
	$(CC) $(FLAGS) -c util.cpp

generator.o :
	$(CC) $(FLAGS) -c generator.cpp

clean :
	rm -f $(OBJS)
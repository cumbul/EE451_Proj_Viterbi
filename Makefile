PROG = main
CC = g++
FLAGS = -std=c++11 -Wall -g
OBJS = main.o viterbi.o parallel_viterbi.o util.o

(PROG) : $(OBJS)
	$(CC) $(FLAGS) -o $(PROG) $(OBJS)
	rm -f $(OBJS)

main.o :
	$(CC) $(FLAGS) -c main.cpp

viterbi.o : viterbi.h
	$(CC) $(FLAGS) -c viterbi.cpp

parallel_viterbi.o : parallel_viterbi.h
	$(CC) $(FLAGS) -c parallel_viterbi.cpp

util.o : util.h
	$(CC) $(FLAGS) -c util.cpp

clean :
	rm -f $(OBJS)
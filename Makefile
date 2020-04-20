PROG = main
CC = g++
FLAGS = -Wall -g
OBJS = main.o viterbi.o parallel_viterbi.o

(PROG) : $(OBJS)
	$(CC) $(FLAGS) -o $(PROG) $(OBJS)
	rm -f $(OBJS)

main.o :
	$(CC) $(FLAGS) -c main.cpp

viterbi.o : viterbi.h
	$(CC) $(FLAGS) -c viterbi.cpp

parallel_viterbi.o : parallel_viterbi.h
	$(CC) $(FLAGS) -c parallel_viterbi.cpp
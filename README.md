# EE451_Proj_Viterbi
Type make to compile files. Type "./main" to run the program.

Parallel Viterbi HPC:

To compile : g++ -std=c++11 -O3 -c main.cpp viterbi.cpp parallel_viterbi.cpp util.cpp -fopenmp

TO link: g++ -std=c++11 -O3 -o main main.o viterbi.o parallel_viterbi.o util.o -fopenmp

TO run: srun -n1 -c2 ./main    (c2 means 2 processors)

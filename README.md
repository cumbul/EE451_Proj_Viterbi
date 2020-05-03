# EE451_Proj_Viterbi
Type make to compile files. Type "./main" to run the program.

The code is dedicated for USC HPC.

# LTDP with openMP
To compile : g++ -std=c++11 -O3 -c main.cpp viterbi.cpp parallel_viterbi.cpp util.cpp -fopenmp

TO link: g++ -std=c++11 -O3 -o main main.o viterbi.o parallel_viterbi.o util.o -fopenmp

TO run: srun -n1 -c2 ./main    (c2 means 2 processors)

# LTDP with MPI
First, setup MPI toolchain by typning ‘source /usr/usc/openmpi/default/setup.sh’

To compile: mpic++ -std=c++11 -O3 -c util.cpp mpi.cpp viterbi.cpp

To link: mpic++ -o mpi util.o mpi.o viterbi.o

TO run: srun -n4 ./mpi (-n4 means 4 nodes)

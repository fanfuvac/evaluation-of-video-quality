SRC=./src/
OBJ=./obj/
all: program

program: cudacode.o
	nvcc  -o program -std=c++11 -L/usr/local/cuda/lib64 -lcuda -lpthread  -Xcompiler -fopenmp  $(SRC)*.cpp  $(OBJ)main.cu.o $(OBJ)SSIM.cu.o  $(OBJ)stvssim.cu.o

cudacode.o:
	
	nvcc -ccbin=/usr/bin/gcc-4.8.5 -std=c++11 -g -G -c  -arch=sm_20 $(SRC)main.cu -o $(OBJ)main.cu.o
	nvcc -ccbin=/usr/bin/gcc-4.8.5 -std=c++11 -g -G -c -arch=sm_20 $(SRC)SSIM.cu -o $(OBJ)SSIM.cu.o
	nvcc  -ccbin=/usr/bin/gcc-4.8.5 -std=c++11 -c  -arch=sm_20  $(SRC)stvssim.cu -o $(OBJ)stvssim.cu.o 

clean: rm -rf *o program

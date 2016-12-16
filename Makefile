SRC=./src/
all: program

program: cudacode.o
	nvcc -o program -g -G -L/usr/local/cuda/lib64 -lcuda -lpthread  -Xcompiler -fopenmp $(SRC)*.cpp  $(SRC)stvssim.cu.o $(SRC)main.cu.o $(SRC)SSIM.cu.o

cudacode.o:
	nvcc -g -G -c -arch=sm_20 $(SRC)stvssim.cu -o $(SRC)stvssim.cu.o
	nvcc -g -G -c -arch=sm_20 $(SRC)main.cu -o $(SRC)main.cu.o
	nvcc -g -G -c -arch=sm_20 $(SRC)SSIM.cu -o $(SRC)SSIM.cu.o

clean: rm -rf *o program

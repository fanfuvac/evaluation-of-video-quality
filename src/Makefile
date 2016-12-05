all: program

program: cudacode.o
	nvcc -o program -g -G -L/usr/local/cuda/lib64 -lcuda -lpthread  -Xcompiler -fopenmp *.cpp  stvssim.cu.o main.cu.o SSIM.cu.o

cudacode.o:
	nvcc -g -c -arch=sm_20 stvssim.cu -o stvssim.cu.o
	nvcc -g -c -arch=sm_20 main.cu -o main.cu.o
	nvcc -g -c -arch=sm_20 SSIM.cu -o SSIM.cu.o

clean: rm -rf *o program

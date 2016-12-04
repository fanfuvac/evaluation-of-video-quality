LIBRARIES       =+ $(_SPACE_) "-L$(TOP)/$(_TARGET_DIR_)/lib$(_TARGET_SIZE_)/stubs" "-L$(TOP)/$(_TARGET_DIR_)/lib$(_TARGET_SIZE_)"

all: program

program: cudacode.o
	g++ -o program -L/usr/local/cuda/lib64 -lcuda -lcudart ./src/*.cpp  cudacode.o 

cudacode.o:
	nvcc -c -arch=sm_20 ./src/*.cu 

clean: rm -rf *o program
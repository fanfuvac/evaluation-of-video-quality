#pragma once
#include <iostream>
#include <fstream>
using namespace std;


__global__ void countRectangleKernel(unsigned char * data1, unsigned char * data2, unsigned char *rects1, unsigned char * rects2, double * out, int size, int width);
double ** countCUDA(FILE ** streams, FILE * ref, int files_count, PictureData * frame, string type, double ** results);
double countSSIM(unsigned char * datain1, unsigned char * datain2, unsigned char * dataC1, unsigned char * dataC2, unsigned char * rects1, unsigned char * rects2, int size, int width, double*& results);



#define RECT_SIZE 64 //number of pixels in rectangle
#define RECT_SQRT 8 //size of side of rectangle
#define C1 6.5025 //(K1*L )**2
#define C2 58.5225 //(K2*L )**2
#define SKIP_SIZE 8 //how many pixels should be skipped when going to next window, should be pow of 2((have to be?) //for counting of vectors, needs to be same as SKIP for 3D-SSIM
//FIXME

extern int THREADS;//how many parallel threads will be started for CUDA kernel


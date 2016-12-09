#pragma once
#include "main.h"   
#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#define RECT_SIZE_ARPS 64 //number of pixels in rectangle
#define RECT_SQRT_ARPS 8 //size of side of rectangle
#define ZERO_MVMT 12 //level of SAD for zero movement prejudgement(2 bit per pixel)
#define FRAME_CNT 33 //count of frames being evaluated by SSIM-3D, need to be odd
#define ROOD_SIZE 2 //size of the rood for ARPS
#define RECT_SQRT_3D 11 //size side of rectangle for 3D-SSIM
#define RECT_SIZE_3D 11*11 //size of rectangle for 3D-SSIM
#define FRAME_SKIP 16 //how mane frames will be skipped between 2 frames (every how many frames is stVSSIM being counted)
#define INT_MAX 2147483647

using namespace std;
//vector used for ARPS
/*struct vector {
	int x;
	int y;
*/
//groups all reusable allocated data for CUDA part of calculation 
struct data_CUDA {
	unsigned char * datain1=NULL;
	unsigned char * datain2=NULL;
	unsigned char * cubes1=NULL;
	unsigned char * cubes2=NULL;
	unsigned char * filters=NULL;
	double * tmpRes=NULL;
	int * filter=NULL;
};
double ** countMetricSTVSSIM_CUDA(FILE ** streams, FILE * ref, int files_count, PictureData * frame, double ** results, int *& frames);

//__global__ void SSIM3DKernel(unsigned char * filters, unsigned char * datain1_CUDA, unsigned char * datain2_CUDA);
double countSTVSSIM_CUDA(unsigned char * datain1, unsigned char * datain2, int size, int width,data_CUDA * dataCuda);
__device__  __host__ void shiftData(unsigned char * data, int size);
__device__ void fillCube(unsigned char * datain, int pos, unsigned char *& out, int width, int height);

__device__ double countMu(unsigned char* filter, unsigned char* cube);
__device__ double countDeltaSqr(unsigned char* filter, unsigned char* cube, double mu);
__device__ double countDelta(unsigned char* filter, unsigned char* cube1, unsigned char* cube2, double muX, double muY);
__device__ double countSSIM3D(unsigned char * filter, unsigned char *  cube1, unsigned char *  cube2);
__device__ void generateCube_CUDA(unsigned char ** cube);

__global__ void SSIM3DKernel(unsigned char * datain1, unsigned char * datain2, unsigned char * cubes1,unsigned char * cubes2,unsigned char * filters,double * tmpRes,int * filter, int width, int height);

__device__  __host__  inline int index(const int x, const int y, const int z, const int aa);
__device__  __host__  inline int indexF(const int x, const int y, int size);
__device__  __host__  inline int index(const int x, const int y, const int z);
__device__  __host__  inline int indexFs(const int x, const int y, const int size, const int file_index);

unsigned char * generateFilters(unsigned char *&);


//vector countARPS(unsigned char * block, unsigned char * framePrev, int x, int y, int width, int height, int T);

int countSAD(unsigned char * rect1, unsigned  char * rect2);

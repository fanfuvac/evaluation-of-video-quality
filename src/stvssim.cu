#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>  
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "stvssim.cuh"
#include "stvssim.h"
#include "SSIM.h"
#include <stdlib.h>     /* abs */
#include <map>
#include <omp.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;


void cudaTest(cudaError_t cudaStatus, string descr) {
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%d", int(cudaStatus));
		fprintf(stderr, strcat("cuda failed: ", descr.c_str()));
		exit(-1);
	}

}

//count STVSSIM metric for all files 
double ** countMetricSTVSSIM_CUDA(FILE ** streams, FILE * ref, int files_count, PictureData * frame, double ** results, int *& frames) {
	cudaError_t cudaStatus;
	
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		return NULL;
		//goto Error;
	}
	//results= new double *[files_count];
#define CUDA_RUNNING 1
	int rec;
	
	size_t pitch;

	//unsigned char * ref_dataHost=new unsigned char frame->size* FRAME_CNT;
	unsigned char * ref_data = new unsigned char[FRAME_CNT*frame->size];
	//cudaStatus = cudaMallocPitch((void **)&ref_data, &pitch, frame->size, FRAME_CNT);
	//cudaTest(cudaStatus, "malloc 1");

	//unsigned char * dataHost = NULL;
	unsigned char * data  = new unsigned char [FRAME_CNT*frame->size];
	//cudaStatus = cudaMallocPitch((void **)&data, &pitch, frame->size, FRAME_CNT*files_count);
	//cudaTest(cudaStatus, "malloc 2");


	unsigned char * dataTrash = new unsigned char[frame->size / 2];
	unsigned char * dataTmp = new unsigned char[frame->size];


	for (int i = FRAME_CNT / 2; i < FRAME_CNT; i++) {
		for (int j = 0; j < files_count; j++) {
			readFromFile(dataTmp, frame->size, streams[j]);
			readFromFile(dataTrash, frame->size / 2, streams[j]);//when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate
			memcpy(data + indexFs(i, 0, frame->size, j), dataTmp, frame->size);
			//cudaMemcpy(data + indexFs(i, 0, frame->size,j),dataTmp,frame->size, cudaMemcpyHostToDevice);
		}
		readFromFile(dataTmp, frame->size, ref);
		readFromFile(dataTrash, frame->size / 2, ref);
		memcpy(ref_data + indexF(i, 0, frame->size), dataTmp, frame->size);
		//cudaMemcpy(ref_data + indexF(i, 0, frame->size), dataTmp, frame->size, cudaMemcpyHostToDevice);
	}

	int i = FRAME_SKIP;
	int j = 0;

	for (; i < frame->frame_count - FRAME_SKIP; i += FRAME_SKIP, j++) {
		shiftData(ref_data, frame->size);
		for (int k = 0; k < files_count; k++) {
			shiftData(data + indexFs(0, 0, frame->size, k), frame->size);
		}
		for (int k = FRAME_CNT / 2 + 1; k < FRAME_CNT; k++) {
			for (int l = 0; l < files_count; l++) {
				readFromFile(dataTmp, frame->size, streams[l]);
				readFromFile(dataTrash, frame->size / 2, streams[l]);//when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate
				memcpy(data + indexFs(i, 0, frame->size, j), dataTmp, frame->size);
				//cudaMemcpy(data + indexFs(k, 0, frame->size,l), dataTmp, frame->size, cudaMemcpyHostToDevice);
			}
			readFromFile(dataTmp, frame->size, ref);
			readFromFile(dataTrash, frame->size / 2, ref);//when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate
			memcpy(ref_data + indexF(i, 0, frame->size), dataTmp, frame->size); //FIXME optimalize
																				//cudaMemcpy(ref_data + indexF(k, 0, frame->size), dataTmp, frame->size, cudaMemcpyHostToDevice);
		}
		double resSSIM, res3D;
		for (int l = 0; l < files_count; l++) {
			res3D = countSTVSSIM_CUDA(ref_data, data + indexFs(0, 0, frame->size, l), frame->size, frame->width);
			//resSSIM = countSSIM(ref_data[FRAME_CNT / 2], data[l][FRAME_CNT / 2], frame->size, frame->width);
			results[l][j] = res3D;
			cout << "3D: " << res3D << " SSIM: " << resSSIM << " Total: " << results[l][j] << endl;

			//cout << j << ": " << results[l][j] << endl;
		}
		//cout << results[j] << endl;

	}
	for (int i = 0; i < files_count; i++) {
		frames[i] = j;
	}
	return results;
}

double countSTVSSIM_CUDA(unsigned char * datain1, unsigned char * datain2, int size, int width) {

	unsigned char * out = new unsigned char[RECT_SIZE];
	int T = ROOD_SIZE;


	unsigned char * filters = new unsigned char[RECT_SIZE_3D*FRAME_CNT * 4];// = new unsigned char ****[CHUNK_SIZE]; //generateFilters();
	double * tmpRes;// = new double[size];
	unsigned char **** cube1 = new unsigned char ***[CHUNK_SIZE]; //generateCube();
	unsigned char **** cube2 = new unsigned char ***[CHUNK_SIZE]; //generateCube();

	int k = 0;
	/*for (int i = 0; i < CHUNK_SIZE; i++) {
	filters[i] = generateFilters();
	cube1[i] = generateCube();
	cube2[i] = generateCube();
	out[i] = new unsigned char[RECT_SIZE];
	}*/

	int rectCount = (size / width - RECT_SQRT_3D)*(width - RECT_SQRT_3D) / SKIP_SIZE / SKIP_SIZE / THREADS*THREADS;
	int blocks = rectCount / THREADS;
	if (rectCount < (size / width - RECT_SQRT_3D)*(width - RECT_SQRT_3D) / SKIP_SIZE / SKIP_SIZE) {
		blocks = rectCount / THREADS + 1;
		rectCount = (size / width - RECT_SQRT_3D)*(width - RECT_SQRT_3D) / SKIP_SIZE / SKIP_SIZE;
	}
	vector vct;
	int * filter = new int[rectCount];
	for (int i = 0; i < size / width - RECT_SQRT_3D; i += SKIP_SIZE) {
		for (int j = 0; j < width - RECT_SQRT_3D; j += SKIP_SIZE) {
			k = (i / SKIP_SIZE)*((width - RECT_SQRT_3D) / SKIP_SIZE) + j / SKIP_SIZE;
			getRect(datain1 + indexF(FRAME_CNT / 2, 0, size), i*width + j, width, out); //FIXME - was i??
																						//if (abs(vct.x) > abs(vct.y)) T = abs(vct.x); FIXME
																						//if (abs(vct.x) < abs(vct.y)) T = abs(vct.y);
			vct = countARPS(out, datain1 + indexF(FRAME_CNT / 2 - 1, 0, size), j, i, width, size / width, T);

			if ((vct.x > vct.y * 2 && vct.x*-1 < 2 * vct.y) || (vct.x < vct.y * 2 && vct.x*-1 > 2 * vct.y)) { //y=0
				filter[k] = 0;
			}
			else if ((vct.y > vct.x * 2 && vct.y > -2 * vct.x) || (vct.y < vct.x * 2 && vct.y < -2 * vct.x)) { //x=0
				filter[k] = 2;
			}
			else if ((vct.y > vct.x / 2 && vct.y < 2 * vct.x) || (vct.y < vct.x / 2 && vct.y > 2 * vct.x)) { //y=x
				filter[k] = 3;
			}
			else if ((vct.y > vct.x / -2 && vct.y < -2 * vct.x) || (vct.y <vct.x / -2 && vct.y>-2 * vct.x)) { //y=-x
				filter[k] = 1;
			}
			else if (vct.x == 0 && vct.y == 0) {
				filter[k] = 8;
			}
			else if (vct.x == vct.y * 2) { //exactly between 2 axes
				filter[k] = 4;
			}
			else if (vct.x == vct.y * -2) {
				filter[k] = 5;
			}
			else if (-2 * vct.x == vct.y) {
				filter[k] = 6;
			}
			else if (2 * vct.x == vct.y) {
				filter[k] = 7;
			}
			else {
				cout << "WUT - nonsense vector " << vct.x << " " << vct.y << endl;
			}

		}
	}

	/*for (int i = 0; i < size / width - RECT_SQRT_3D; i += SKIP_SIZE) {
	for (int j = 0; j < width - RECT_SQRT_3D; j += SKIP_SIZE) {
	generateCube_CUDA()
	fillCube(datain1, i*width + j, cubeTmp, width);
	}
	}*/

	cudaError_t cudaStatus;
	size_t pitch;
	unsigned char * filters_CUDA=NULL;
	int * filter_CUDA=NULL;
	unsigned char * datain1_CUDA=NULL;
	unsigned char * datain2_CUDA=NULL;
	cudaStatus = cudaMallocPitch((void **)&filters_CUDA, &pitch, RECT_SIZE_3D, FRAME_CNT);
	cudaStatus = cudaMallocPitch((void **)&datain1_CUDA, &pitch, width, size / width);
	cudaStatus = cudaMallocPitch((void **)&datain2_CUDA, &pitch, width, size / width);
	cudaStatus = cudaMallocPitch((void **)&tmpRes, &pitch, sizeof(double), rectCount);
	cudaTest(cudaStatus, "malloc filters");
	generateFilters(filters);
	cudaStatus = cudaMemcpy((void*)&filters_CUDA, (const void*)filters, FRAME_CNT*RECT_SIZE_3D, cudaMemcpyHostToDevice);
cudaTest(cudaStatus, "memcpy filters");
	cudaStatus = cudaMemcpy((void*)&datain1_CUDA, (void*)datain1, size, cudaMemcpyHostToDevice);
	cudaTest(cudaStatus, "memcpy data1");
	cudaStatus = cudaMemcpy((void*)&datain2_CUDA, (void*)datain2, size, cudaMemcpyHostToDevice);
	cudaTest(cudaStatus, "memcpy data2");
	cudaStatus = cudaMemcpy((void*)&filter_CUDA, (void*)filter, rectCount * sizeof(int), cudaMemcpyHostToDevice);
	cudaTest(cudaStatus, "memcpy filter");

	SSIM3DKernel << <blocks, THREADS >> > (filters_CUDA, datain1_CUDA, datain2_CUDA, filter_CUDA, tmpRes, width, size / width);

	k = (size / width - RECT_SQRT) / SKIP_SIZE*(width - RECT_SQRT) / SKIP_SIZE;
	double * tmpRes2 = new double[k];
	cudaMemcpy((void*)tmpRes2, (void*)tmpRes2, k*sizeof(double), cudaMemcpyDeviceToHost);
	
	
	double res = countRes(tmpRes2, k);
	delete[] tmpRes2;

	/*
	for (int l = 0; l < CHUNK_SIZE; l++) {
	for (int i = 0; i < FRAME_CNT; i++) {
	for (int j = 0; j < RECT_SQRT_3D; j++) {

	delete[] cube1[l][i][j];
	delete[] cube2[l][i][j];
	//delete[] filters[l][0][j];
	}
	delete[] cube1[l][i];
	delete[] cube2[l][i];
	//delete[] filters[l][i];
	}
	delete[] cube1[l];
	delete[] cube2[l];
	//delete[] filters[l];
	delete[] out[l];
	}
	delete[] cube1;
	delete[] cube2;
	//delete[] filters;
	delete[] out;

	for (int l = 0; l < CHUNK_SIZE; l++) {
	for (int i = 0; i < 4; i++) {
	for (int j = 0; j < FRAME_CNT; j++) {
	for (int k = 0; k < RECT_SQRT_3D; k++) {
	delete[] filters[l][i][j][k];
	}
	delete[] filters[l][i][j];
	}
	delete[] filters[l][i];
	}
	delete[] filters[l];
	}
	delete[] filters;
	*/

	return res;

}
__global__ void SSIM3DKernel(unsigned char * filters, unsigned char * datain1_CUDA, unsigned char * datain2_CUDA, int * filter, double * results, int width, int height) {
	int i = threadIdx.x;
	int j = blockIdx.x;
	int pos = j*THREADS + i;
	unsigned char *cube1;
	cube1 = generateCube_CUDA(cube1);
	unsigned char *cube2;
	cube2 = generateCube_CUDA(cube2);
	fillCube(datain1_CUDA, j*width + i, cube1, width, height);
	fillCube(datain2_CUDA, j*width + i, cube2, width, height);

	double res0 = countSSIM3D(filters, cube1, cube2);
	double res1 = countSSIM3D(filters + FRAME_CNT*RECT_SIZE_3D * 1, cube1, cube2);
	double res2 = countSSIM3D(filters + FRAME_CNT*RECT_SIZE_3D * 2, cube1, cube2);
	double res3 = countSSIM3D(filters + FRAME_CNT*RECT_SIZE_3D * 3, cube1, cube2);
	switch (filter[pos]) {
	case 0:
		results[pos] = res0;
	case 1:
		results[pos] = res1;
	case 2:
		results[pos] = res2;
	case 3:
		results[pos] = res3;
	case 4:
		results[pos] = (res0 + res3) / 2;
	case 5:
		results[pos] = (res0 + res1) / 2;
	case 6:
		results[pos] = (res1 + res2) / 2;
	case 7:
		results[pos] = (res2 + res3) / 2;
	case 8:
		results[pos] = (res0 + res1 + res2 + res3) / 4;
	}


}

__device__ double countSSIM3D(unsigned char * filter, unsigned char *  cube1, unsigned char *  cube2) {
	double muX = countMu(filter, cube1);
	double muY = countMu(filter, cube2);

	double deltaSqrX = countDeltaSqr(filter, cube1, muX);
	double deltaSqrY = countDeltaSqr(filter, cube2, muY);

	double delta = countDelta(filter, cube1, cube2, muX, muY);

	double ssim3D = ((2 * muX*muY + C1)*(2 * delta + C2)) / ((muX*muX + muY*muY + C1)*(deltaSqrX + deltaSqrY + C2));
	return ssim3D;
}
__device__ double countMu(unsigned char* filter, unsigned char* cube) {
	double res = 0;
	for (int alpha = 0; alpha < RECT_SQRT_3D; alpha++) {
		for (int beta = 0; beta < RECT_SQRT_3D; beta++) {
			for (int gamma = 0; gamma < FRAME_CNT; gamma++) {

				res += filter[index(gamma, alpha, beta)] * cube[index(gamma, alpha, beta)];
			}
		}
	}
	return res / (RECT_SQRT_3D*FRAME_CNT);
}

__device__ double countDeltaSqr(unsigned char* filter, unsigned char* cube, double mu) {
	double res = 0;
	for (int alpha = 0; alpha < RECT_SQRT_3D; alpha++) {
		for (int beta = 0; beta < RECT_SQRT_3D; beta++) {
			for (int gamma = 0; gamma < FRAME_CNT; gamma++) {
				res += filter[index(gamma, alpha, beta)] * (cube[index(gamma, alpha, beta)] - mu)*(cube[index(gamma, alpha, beta)] - mu);
			}
		}
	}
	return res / (RECT_SQRT_3D*FRAME_CNT);
}


__device__ double countDelta(unsigned char* filter, unsigned char* cube1, unsigned char* cube2, double muX, double muY) {
	double res = 0;
	for (int alpha = 0; alpha < RECT_SQRT_3D; alpha++) {
		for (int beta = 0; beta < RECT_SQRT_3D; beta++) {
			for (int gamma = 0; gamma < FRAME_CNT; gamma++) {
				res += filter[index(gamma, alpha, beta)] * (cube1[index(gamma, alpha, beta)] - muX)*(cube2[index(gamma, alpha, beta)] - muX);
			}
		}
	}
	return res / (RECT_SQRT_3D*FRAME_CNT);
}

//Generates 3D array used for SSIM 3D
__device__ unsigned char * generateCube_CUDA(unsigned char* cube) {
	cudaError_t cudaStatus;
	//	unsigned char* cube;
	size_t pitch;
	cube = new unsigned char[FRAME_CNT*RECT_SIZE_3D];
	/*cudaStatus = cudaMallocPitch((void **)&cube, &pitch, FRAME_CNT, RECT_SIZE_3D);
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMalloc failed 3!");
	return NULL;
	}*/

	return cube;
}

//Fill 3D array with data of surroundings pixels
__device__ void fillCube(unsigned char * datain, int pos, unsigned char * out, int width, int height) {
	for (int i = 0; i < FRAME_CNT; i++) {
		for (int j = 0; j < RECT_SQRT_3D; j++) {
			for (int k = 0; k < RECT_SQRT_3D; k++) {
				out[index(i, j, k)] = datain[indexF(i, pos + width*j, width*height)];
				//memcpy(out[i][j], datain[i] + pos + width*j, RECT_SQRT_3D);
			}
		}
	}
}

//generate cube filters, vertical, horizontal and 2 inclined are being created
//receives already malloced referenco to an CUDA array
__host__ unsigned char * generateFilters(unsigned char*& out) {

	for (int i = 0; i < FRAME_CNT; i++) { //FIXME modify to create 2D array and memcpy
		for (int j = 0; j < RECT_SQRT_3D; j++) {
			for (int k = 0; k < RECT_SQRT_3D; k++) {
				if (j == RECT_SQRT_3D / 2) { //horizontal filter
					out[index(0, i, j, k)] = 1;
				}
				else {
					out[index(0, i, j, k)] = 0;
				}
				if (k == RECT_SQRT_3D / 2) { //vertical filter
					out[index(1, i, j, k)] = 1;
				}
				else {
					out[index(1, i, j, k)] = 0;
				}


				if (j == k) { //x=y filter
					out[index(2, i, j, k)] = 1;
				}
				else {
					out[index(2, i, j, k)] = 0;
				}
				if (k + j == RECT_SQRT_3D) { //x=-y filter
					out[index(3, i, j, k)] = 1;
				}
				else {
					out[index(3, i, j, k)] = 0;
				}
			}
		}
	}
	return out;
}



/*
vector countARPS(unsigned char * block, unsigned char * framePrev, int x, int y, int width, int height, int T) {
unsigned char * out = new unsigned char[RECT_SIZE];
getRect(framePrev, x*y, width, out);
int sad = countSAD(block, out);
vector vOut;
if (sad < ZERO_MVMT) {
vOut.x = 0;
vOut.y = 0;
return vOut;
}
map<int, int > past;
int xOrig = x;
int yOrig = y;
int res[5];
while (1) {
getRect(framePrev, x + y*width, width, out);
res[0] = countSAD(block, out);
if (x - T - RECT_SQRT / 2 > 0) {
getRect(framePrev, (x - T) + y*width, width, out);
res[3] = countSAD(block, out);
}
else res[3] = INT_MAX;
if (x + T + RECT_SQRT / 2< width) {
getRect(framePrev, (x + T) + y*width, width, out);
res[1] = countSAD(block, out);
}
else res[1] = INT_MAX;
if (y + T + RECT_SQRT / 2< height) {
getRect(framePrev, x + (y + T)*width, width, out);
res[2] = countSAD(block, out);
}
else res[2] = INT_MAX;
if (y - T - RECT_SQRT / 2 > 0) {
getRect(framePrev, x + (y - T)*width, width, out);
res[4] = countSAD(block, out);
}
else res[4] = INT_MAX;
int min = INT_MAX;
int minPos = 0;
for (int i = 0; i < 5; i++) {
if (res[i] < min) {
min = res[i];
minPos = i;
}
}
if (minPos == 0) {
vOut.x = x - xOrig;
vOut.y = y - yOrig;
delete[] out;
return vOut;
}
switch (minPos) {
case 1:
x += T;
break;
case 2:
y += T;
break;
case 3:
x -= T;
break;
case 4:
y -= T;
break;
}
}
}
*/
__device__ __host__ void shiftData(unsigned char * data, int size) {
	for (int i = 0; i < FRAME_CNT / 2 + 1; i++) {
#if defined(CUDA_RUNNING)
		cudaMemcpy(data + indexF(i, 0, size), data + indexF(i + FRAME_CNT / 2, 0, size), size, cudaMemcpyDeviceToDevice);
#else
		memcpy(data + indexF(i, 0, size), data + indexF(i + FRAME_CNT / 2, 0, size), size);
#endif
	}

}

__device__ __host__ int countSAD(unsigned char * rect1, unsigned  char * rect2) {
	int sad = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		sad += abs(rect1[i] - rect2[i]);
	}
	return sad;
}

//used for flattening of arrays because CUDA cannot work with multi dimensional arrays, FRAME_CNT frames
__device__  __host__  inline int indexF(const int x, const int y, const int size) {
	return x * size + y;
}

//used for flattening of arrays because CUDA cannot work with multi dimensional arrays, FRAME_CNT frames with files_count files
__device__  __host__  inline int indexFs(const int x, const int y, const int size, const int file_index) {
	return x * size *file_index + y;
}

//used for flattening of arrays because CUDA cannot work with multi dimensional arrays, cube
__device__  __host__ inline int index(const int x, const int y, const int z) {
	return x * RECT_SIZE_3D + y * RECT_SQRT_3D + z;
}

//used for flattening of arrays because CUDA cannot work with multi dimensional arrays, filter
__device__  __host__ inline int index(const int x, const int y, const int z, const int aa) {
	return x * RECT_SIZE_3D*FRAME_CNT + y * RECT_SIZE_3D + z *RECT_SQRT_3D + aa;
}

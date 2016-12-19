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
#include "ssim.cuh"
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

///count STVSSIM metric for all files 
double ** countMetricSTVSSIM_CUDA(FILE ** streams, FILE * ref, int files_count, PictureData * frame, double ** results, int *& frames) {
	cudaError_t cudaStatus;

	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		exit(-1);
		return NULL;
	}
	int rec;
	size_t pitch;

	unsigned char * ref_data = new unsigned char[FRAME_CNT*frame->size];
	unsigned char * data = new unsigned char[FRAME_CNT*frame->size*files_count];
	unsigned char * dataTrash = new unsigned char[frame->size / 2];
	unsigned char * dataTmp = new unsigned char[frame->size];

	for (int i = FRAME_CNT / 2; i < FRAME_CNT; i++) {
		for (int j = 0; j < files_count; j++) {
			readFromFile(dataTmp, frame->size, streams[j]);
			readFromFile(dataTrash, frame->size / 2, streams[j]);//when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate
			memcpy(data + indexFs(i, 0, frame->size, j), dataTmp, frame->size);
		}
		readFromFile(dataTmp, frame->size, ref);
		readFromFile(dataTrash, frame->size / 2, ref);
		memcpy(ref_data + indexF(i, 0, frame->size), dataTmp, frame->size);
	}

	int i = FRAME_SKIP;
	int j = 0;
	data_CUDA * dataCuda = new data_CUDA;

	cudaStatus = cudaMallocPitch((void **)&dataCuda->filters, &pitch, RECT_SIZE_3D, FRAME_CNT * 4);//4 different filters with size of one cube
	cudaTest(cudaStatus, "malloc filters");
	cudaStatus = cudaMallocPitch((void **)&dataCuda->datain1, &pitch, FRAME_CNT, frame->size);
	cudaTest(cudaStatus, "malloc datain1");
	cudaStatus = cudaMallocPitch((void **)&dataCuda->datain2, &pitch, FRAME_CNT, frame->size);
	int k = (frame->size / frame->width - RECT_SQRT) / SKIP_SIZE*(frame->width - RECT_SQRT) / SKIP_SIZE;
	cudaStatus = cudaMallocPitch((void **)&dataCuda->tmpRes, &pitch, sizeof(double), k);
	cudaTest(cudaStatus, "malloc datain2");
	cudaStatus = cudaMallocPitch((void **)&dataCuda->filter, &pitch, sizeof(int), k);
	cudaTest(cudaStatus, "malloc filter");

	cudaStatus = cudaMallocPitch((void **)&dataCuda->cubes1, &pitch, RECT_SIZE_3D*FRAME_CNT, k);
	cudaTest(cudaStatus, "malloc filter");
	cudaStatus = cudaMallocPitch((void **)&dataCuda->cubes2, &pitch, RECT_SIZE_3D*FRAME_CNT, k);
	cudaTest(cudaStatus, "malloc filter");

	//preparation for SSIM GPU part
	unsigned char * dataRef;
	cudaStatus = cudaMalloc((void **)&dataRef, frame->width*frame->height);
	cudaTest(cudaStatus, "malloc dataRef");

	unsigned char ** dataGPU = new unsigned char *[files_count];
	//allocated the device memory for source array  
	for (int i = 0; i < files_count; i++) {
		cudaStatus = cudaMalloc((void **)&dataGPU[i], frame->width*frame->height);
		cudaTest(cudaStatus, "malloc data GPU SSIM");
	}

	unsigned char ** rects = new unsigned char *[files_count];
	for (int i = 0; i < 2; i++) { //if replaced with files_count it can start kernels in parallel
		cudaStatus = cudaMallocPitch((void **)&rects[i], &pitch, RECT_SIZE, frame->width*frame->height / SKIP_SIZE / SKIP_SIZE);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed 3!");
			return NULL;
		}
	}

	double * resultsFrame;
	cudaStatus = cudaMallocPitch((void**)&resultsFrame, &pitch, sizeof(double), frame->width*frame->height / SKIP_SIZE / SKIP_SIZE);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed 5!");

		return NULL;
	}
	//end of preparation for GPU SSIM

	for (; i < frame->frame_count - FRAME_SKIP-1; i += FRAME_SKIP, j++) {
		shiftData(ref_data, frame->size);
		for (int k = 0; k < files_count; k++) {
			shiftData(data + indexFs(0, 0, frame->size, k), frame->size);
		}
		for (int k = FRAME_CNT / 2 + 1; k < FRAME_CNT; k++) {
			for (int l = 0; l < files_count; l++) {
				readFromFile(dataTmp, frame->size, streams[l]);
				readFromFile(dataTrash, frame->size / 2, streams[l]);//when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate
				memcpy(data + indexFs(k, 0, frame->size, l), dataTmp, frame->size);
				//cudaMemcpy(data + indexFs(k, 0, frame->size,l), dataTmp, frame->size, cudaMemcpyHostToDevice);
			}
			readFromFile(dataTmp, frame->size, ref);
			readFromFile(dataTrash, frame->size / 2, ref);//when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate
			memcpy(ref_data + indexF(k, 0, frame->size), dataTmp, frame->size); //FIXME optimalize
																				//cudaMemcpy(ref_data + indexF(k, 0, frame->size), dataTmp, frame->size, cudaMemcpyHostToDevice);
		}
		double resSSIM=0, res3D;
		cout << "file" << "\t" << "3D" << "\t" << "SSIM" << "\t" << "total" << endl;
		for (int l = 0; l < files_count; l++) {
			res3D = countSTVSSIM_CUDA(ref_data, data + indexFs(0, 0, frame->size, l), frame->size, frame->width, dataCuda);
			
			//now we can either use CPU or GPU version of SSIM algorithm
			//resSSIM = countSSIM(ref_data+indexF(FRAME_CNT / 2,0,frame->size), data+ indexFs(FRAME_CNT / 2, 0, frame->size, l), frame->size, frame->width);
			resSSIM=countSSIM(ref_data + indexF(FRAME_CNT / 2, 0, frame->size), data + indexFs(FRAME_CNT / 2, 0, frame->size, l), dataRef, dataGPU[l], rects[0], rects[1], frame->size, frame->width, resultsFrame);
			results[l][j] = res3D*resSSIM;
			cout << l << "\t" << res3D << "\t" << resSSIM << "\t" << results[l][j] << endl;
			//cout << "3D: " << res3D << " SSIM: " << resSSIM << " Total: " << results[l][j] << endl;
		}
	}
	for (int i = 0; i < files_count; i++) {
		frames[i] = j;
	}
	return results;
}

/*!counting STVSSIM for given set of frames
 @param[in] datain1,datain2 input datasets
 @param[in] size size in bytes of one image
 @param[in] width width of one image
 @param[in] dataCuda struct with preallocated space on CUDA device
*/
double countSTVSSIM_CUDA(unsigned char * datain1, unsigned char * datain2, int size, int width, data_CUDA * dataCuda) {

	unsigned char * out = new unsigned char[RECT_SIZE];
	int T = ROOD_SIZE;


	unsigned char * filters = new unsigned char[RECT_SIZE_3D*FRAME_CNT * 4];// = new unsigned char ****[CHUNK_SIZE]; //generateFilters();

	//unsigned char **** cube1 = new unsigned char ***[CHUNK_SIZE]; //generateCube();
	//unsigned char **** cube2 = new unsigned char ***[CHUNK_SIZE]; //generateCube();

	int k = 0;

	int rectCount = (int)((size / width - RECT_SQRT_3D) / SKIP_SIZE + 1)*	(int)((width - RECT_SQRT_3D) / SKIP_SIZE + 1) / THREADS*THREADS;
	int blocks = rectCount / THREADS;
	if (rectCount < (int)((size / width - RECT_SQRT_3D) / SKIP_SIZE + 1)*	(int)((width - RECT_SQRT_3D) / SKIP_SIZE + 1)) {
		blocks = rectCount / THREADS + 1;
		rectCount = (int)((size / width - RECT_SQRT_3D) / SKIP_SIZE + 1)*	(int)((width - RECT_SQRT_3D) / SKIP_SIZE + 1);
	}
	vector vct;
	int * filter = new int[rectCount];
	for (int i = 0; i < size / width - RECT_SQRT_3D; i += SKIP_SIZE) {
		for (int j = 0; j < width - RECT_SQRT_3D; j += SKIP_SIZE, k++) {
			getRect(datain1 + indexF(FRAME_CNT / 2, 0, size), i*width + j, width, out); 
			
			vct = countARPS(out, datain1 + indexF(FRAME_CNT / 2 - 1, 0, size), j, i, width, size / width, T);
			//if (abs(vct.x) > abs(vct.y)) T = abs(vct.x); FIXME
			//if (abs(vct.x) < abs(vct.y)) T = abs(vct.y);
			
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
			//printf("vct: %d %d %d %d\n",vct.x,vct.y,k,filter[k]);
		}
	}

	cudaError_t cudaStatus;
	size_t pitch;
	
	generateFilters(filters);

	cudaStatus = cudaMemcpy((void*)dataCuda->filters, (const void*)filters, FRAME_CNT*RECT_SIZE_3D * 4, cudaMemcpyHostToDevice);
	cudaTest(cudaStatus, "memcpy filters");
	cudaStatus = cudaMemcpy((void*)dataCuda->datain1, (void*)datain1, size*FRAME_CNT, cudaMemcpyHostToDevice);
	cudaTest(cudaStatus, "memcpy data1");
	cudaStatus = cudaMemcpy((void*)dataCuda->datain2, (void*)datain2, size*FRAME_CNT, cudaMemcpyHostToDevice);
	cudaTest(cudaStatus, "memcpy data2");
	cudaStatus = cudaMemcpy((void*)dataCuda->filter, (void*)filter, rectCount * sizeof(int), cudaMemcpyHostToDevice);
	cudaTest(cudaStatus, "memcpy filter");


	//dont know how to transfer struct to device
	SSIM3DKernel << <blocks, THREADS >> > (dataCuda->datain1, dataCuda->datain2, dataCuda->cubes1, dataCuda->cubes2, dataCuda->filters, dataCuda->tmpRes, dataCuda->filter, width, size / width);

	cudaDeviceSynchronize();
	double * tmpRes2 = new double[rectCount];
	cudaMemcpy(tmpRes2, dataCuda->tmpRes, rectCount * sizeof(double), cudaMemcpyDeviceToHost);


	double res = countRes(tmpRes2, rectCount);
	delete[] tmpRes2;
	delete[] filters;

	return res;

}

/*!Kernel which is being executed on CUDA device, sounting SSIM-3D part
 @param[in] datain1,datain2 input datasets, allocated on GPU card
 @param cubes1,cubes2 allocated space for data, orriginaly as 3D array, here only 1D array due to CUDA
 @param filters similar to cubes1, used for filtering the data
 @param tmpRes preallocated space on CUDA device used for temporary stored results
 @param filter list of filters being used for each rectangle
 @param width,height size of respective dimension of the frame
*/
__global__ void SSIM3DKernel(unsigned char * datain1, unsigned char * datain2, unsigned char * cubes1, unsigned char * cubes2, unsigned char * filters, double * tmpRes, int * filter, int width, int height) {
	int i = threadIdx.x;
	int j = blockIdx.x;
	int pos = j*blockDim.x + i;
	if (pos >= ((int)((height - RECT_SQRT_3D) / SKIP_SIZE + 1)*	(int)((width - RECT_SQRT_3D) / SKIP_SIZE + 1))) {
		return;
	}
	unsigned char * cube1 = cubes1 + RECT_SIZE_3D*FRAME_CNT*pos;
	unsigned char *cube2 = cubes2 + RECT_SIZE_3D*FRAME_CNT*pos;
	int leftover = (width - RECT_SQRT_3D) % SKIP_SIZE;
	int line = width - RECT_SQRT_3D - leftover+ SKIP_SIZE; //effective size of the line
	int start = (pos*SKIP_SIZE) % ((line)) + (pos*SKIP_SIZE) / line * width*SKIP_SIZE;
	//printf("%d %d: %d ",blockIdx.x, threadIdx.x, start);
	fillCube(datain1, start, cube1, width, height);
	fillCube(datain2, start, cube2, width, height);

	double res0 = countSSIM3D(filters, cube1, cube2);
	double res1 = countSSIM3D(filters + FRAME_CNT*RECT_SIZE_3D * 1, cube1, cube2);
	double res2 = countSSIM3D(filters + FRAME_CNT*RECT_SIZE_3D * 2, cube1, cube2);
	double res3 = countSSIM3D(filters + FRAME_CNT*RECT_SIZE_3D * 3, cube1, cube2);
	//printf("%d %d %d: %f %f %f %f %d\n", start, blockIdx.x, threadIdx.x, res0, res1, res2, res3, filter[pos]);
	switch (filter[pos]) {
	case 0:
		tmpRes[pos] = res0;
		break;
	case 1:
		tmpRes[pos] = res1;
		break;
	case 2:
		tmpRes[pos] = res2;
		break;
	case 3:
		tmpRes[pos] = res3;
		break;
	case 4:
		tmpRes[pos] = (res0 + res3) / 2;
		break;
	case 5:
		tmpRes[pos] = (res0 + res1) / 2;
		break;
	case 6:
		tmpRes[pos] = (res1 + res2) / 2;
		break;
	case 7:
		tmpRes[pos] = (res2 + res3) / 2;
		break;
	case 8:
		tmpRes[pos] = (res0 + res1 + res2 + res3) / 4;
		break;
	}


}


/*!Computation of SSIM-3D part, computes one cube, being executed on device
 @param[in] filter respective filter which will be used for this function
 @param cubes1,cubes2 data with cumputed part of framesorriginaly as 3D array, here only 1D array due to CUDA
*/
__device__ double countSSIM3D(unsigned char * filter, unsigned char *  cube1, unsigned char *  cube2) {
	double muX = countMu(filter, cube1);
	//printf("%d %d: %.0f ",blockIdx.x, threadIdx.x, muX);
	double muY = countMu(filter, cube2);

	double deltaSqrX = countDeltaSqr(filter, cube1, muX);
	double deltaSqrY = countDeltaSqr(filter, cube2, muY);

	double delta = countDelta(filter, cube1, cube2, muX, muY);

	double ssim3D = ((2 * muX*muY + C1)*(2 * delta + C2)) / ((muX*muX + muY*muY + C1)*(deltaSqrX + deltaSqrY + C2));
	return ssim3D;
}
__device__ double countMu(unsigned char* filter, unsigned char* cube) {
	double res = 0, res2 = 0;
	for (int alpha = 0; alpha < RECT_SQRT_3D; alpha++) {
		for (int beta = 0; beta < RECT_SQRT_3D; beta++) {
			for (int gamma = 0; gamma < FRAME_CNT; gamma++) {
				res2 += filter[index(gamma, alpha, beta)];
				res += filter[index(gamma, alpha, beta)] * cube[index(gamma, alpha, beta)];
				//printf("%.0f ",res);

			}
		}
		//printf("res: %d\n",res);
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

///Generates 3D array used for SSIM 3D
__device__ void generateCube_CUDA(unsigned char** cube) {
	cudaError_t cudaStatus;
	//	unsigned char* cube;
	size_t pitch;
	*cube = (unsigned char*)malloc(FRAME_CNT*RECT_SIZE_3D);


	//return cube;
}

///Fill 3D array with data of surroundings pixels
__device__ void fillCube(unsigned char * datain, int pos, unsigned char *& out, int width, int height) {
	for (int i = 0; i < FRAME_CNT; i++) {
		for (int j = 0; j < RECT_SQRT_3D; j++) {
			for (int k = 0; k < RECT_SQRT_3D; k++) {
				out[index(i, j, k)] = datain[indexF(i, pos + width*j + k, width*height)];
			}
		}
	}
}
/*!generate cube filters, vertical, horizontal and 2 inclined are being created
 @param out already malloced referenco to an CUDA array
*/
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
//implement this for CUDA? Probably would not be effective because of too many if and splitting of the code
vector countARPS(unsigned char * block, unsigned char * framePrev, int x, int y, int width, int height, int T) {
}
*/
__device__ __host__ void shiftData(unsigned char * data, int size) {
	for (int i = 0; i < FRAME_CNT / 2 + 1; i++) {
#if defined(__CUDA_ARCH__)
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

///used for flattening of arrays because CUDA cannot work with multi dimensional arrays, FRAME_CNT frames
__device__  __host__  inline int indexF(const int x, const int y, const int size) {
	return x * size + y;
}

///used for flattening of arrays because CUDA cannot work with multi dimensional arrays, FRAME_CNT frames with files_count files
__device__  __host__  inline int indexFs(const int x, const int y, const int size, const int file_index) {
	return file_index * size * FRAME_CNT + x * size + y;
}

///used for flattening of arrays because CUDA cannot work with multi dimensional arrays, cube
__device__  __host__ inline int index(const int x, const int y, const int z) {
	return x * RECT_SIZE_3D + y * RECT_SQRT_3D + z;
}

///used for flattening of arrays because CUDA cannot work with multi dimensional arrays, filter
__device__  __host__ inline int index(const int x, const int y, const int z, const int aa) {
	return x * RECT_SIZE_3D*FRAME_CNT + y * RECT_SIZE_3D + z *RECT_SQRT_3D + aa;
}

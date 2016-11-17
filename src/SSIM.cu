#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>  
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "main.h"  
#include "ssim.cuh"  
//#include "..\..\Video_comparsion\Video_comparsion\PSNR.h"   
#include <omp.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <errno.h>
#include <stdio.h>
using namespace std;
//cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);


__host__ __device__ double countAvg(unsigned char * data);
__host__ __device__ double countVariance(unsigned char * data, double avg);
__host__ __device__ double countCovariance(unsigned char * data1, unsigned char * data2, double avg1, double avg2);
__host__ __device__ void getRect(unsigned char* data, int start, int width, unsigned char * out);
__host__ __device__ double countRectangle(unsigned char * data1, unsigned char * data2);
double countRes(double * tmpRes, int count);

double countSSIM(unsigned char * datain1, unsigned char * datain2, unsigned char * dataC1, unsigned char * dataC2, unsigned char * rects1, unsigned char * rects2, int size, int width, double*& results) {
	//unsigned char * data1 = (unsigned char*)datain1;
	//unsigned char * data2 = (unsigned char*)datain2;
	cudaError_t cudaStatus;
	//double * tmpRes = new double[size];


	/*
	if (data1==0 or data2==0){
	//return NULL;
	//cout<<"error in allocation"<<endl;
	}*/
	/*getLuma(datain1, data1, size);
	getLuma(datain2, data2, size);*/
	/*unsigned char * rect1 = new unsigned char[RECT_SIZE];
	unsigned char * rect2 = new unsigned char[RECT_SIZE];
	int k = 0;*/



	//#pragma omp parallel
	//nthreads = omp_get_num_threads();


	cudaStatus = cudaMemcpy(dataC1, datain1, size, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		printf("%s\n", cudaGetErrorString(cudaStatus));
		return NULL;
	}
	cudaStatus = cudaMemcpy(dataC2, datain2, size, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		printf("%s\n", cudaGetErrorString(cudaStatus));
		return NULL;
	}


	//	#pragma omp parallel for schedule(static, 20)
	int rectCount = (size / width - RECT_SQRT)*(width - RECT_SQRT) / SKIP_SIZE / SKIP_SIZE / THREADS*THREADS;
	int blocks = rectCount / THREADS;
	if (rectCount< (size / width - RECT_SQRT)*(width - RECT_SQRT) / SKIP_SIZE / SKIP_SIZE) {
		blocks = rectCount / THREADS + 1;
		rectCount = (size / width - RECT_SQRT)*(width - RECT_SQRT) / SKIP_SIZE / SKIP_SIZE;
	}
	countRectangleKernel << <blocks, THREADS >> >(dataC1, dataC2, rects1, rects2, results, size, width); //FIXME - need to adjust size to count up to THREADS last rectangles!!
																										 //cudaDeviceSynchronize();
																										 /*for (int i = 0; i < size / width - RECT_SQRT; i += SKIP_SIZE) {

																										 for (int j = 0; j < width - RECT_SQRT; j += SKIP_SIZE, k++) {
																										 //for (int i = 0; i < size-(RECT_SQRT-1)*width; i+=SKIP_SIZE) {

																										 //if (tmpRes[k] < 0) cout << "low result: " << i<< ": " << j<< " :" << tmpRes[k] << endl;

																										 }
																										 }


																										 double res = countRes(tmpRes, k);
																										 */



	double * resultsOut = new double[rectCount];
	cudaStatus = cudaMemcpy(resultsOut, results, rectCount * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		printf("%s\n", cudaGetErrorString(cudaStatus));
		return NULL;
	}
	double output = countRes(resultsOut, rectCount);

	delete resultsOut;
	//cout << output << endl;
	return output;
}





//return one rectangle with RECT_SIZE pixels
__host__ __device__ void getRect(unsigned char* data, int start, int width, unsigned char * out) {
#if defined(__CUDA_ARCH__)
	for (int i = 0; i<RECT_SQRT; i++) {
		for (int j = 0; j < RECT_SQRT; j++) {
			out[i*RECT_SQRT + j] = data[start + i*width + j];
		}
		//cudaMemcpy(out + i*RECT_SQRT, data + start + i*width, RECT_SQRT, cudaMemcpyDeviceToDevice);
	}
#else
	for (int i = 0; i<RECT_SQRT; i++) {
		memcpy(out + i*RECT_SQRT, data + start + i*width, RECT_SQRT);
	}
#endif
	//return out;
}

//count ssim of one rectangle with RECT_SIZE pixels
__host__ __device__ double countRectangle(unsigned char * data1, unsigned char * data2) {

	double avg1 = countAvg(data1);
	double avg2 = countAvg(data2);

	double var1 = countVariance(data1, avg1);
	double var2 = countVariance(data2, avg2);

	double cov = countCovariance(data1, data2, avg1, avg2);


	double ssim = ((2 * avg1*avg2 + C1)*(2 * cov + C2)) / ((avg1*avg1 + avg2*avg2 + C1)*(var1 + var2 + C2));
	return ssim;
}
//count avg value of given rectangle 
__host__  __device__ double countAvg(unsigned char * data) {
	double avg = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		avg += data[i];
	}
	avg = avg / (double)RECT_SIZE;
	return avg;
}

//count variance of given rectangle
__host__ __device__ double countVariance(unsigned char * data, double avg) {
	double var = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		var += (data[i] - avg)*(data[i] - avg);
	}
	var = var / (double)RECT_SIZE;
	return var;
}

//count covariance of given rectangle
__host__ __device__ double countCovariance(unsigned char * data1, unsigned char * data2, double avg1, double avg2) {
	double cov = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		cov += (data1[i] - avg1)*(data2[i] - avg2);
	}
	cov = cov / (double)RECT_SIZE;
	//if (cov < 0) cout << "neg "<<cov << endl;
	return cov;
}




//count average SSIM value from SSIM values per rectangle
double countRes(double * tmpRes, int count) {
	double sum = 0;
	for (int i = 0; i < count; i += 1) {
		//cout << tmpRes[i]<<endl;
		sum += tmpRes[i];

	}
	//cout << "frame" << endl;
	return sum / (double)count;

}
__global__ void countRectangleKernel(unsigned char * data1, unsigned char * data2, unsigned char * rects1, unsigned char * rects2, double * out, int size, int width) {
	int i = threadIdx.x;
	int j = blockIdx.x;
	int a;
	int pos = j*THREADS + i;
	int x = (pos*SKIP_SIZE) % (width - RECT_SQRT);
	int y = (pos*SKIP_SIZE) / (width - RECT_SQRT)*SKIP_SIZE;
	if (pos<size / SKIP_SIZE / SKIP_SIZE) {
		getRect(data1, x+y*width, width, rects1 + (pos)*RECT_SIZE);
		getRect(data2, x + y*width, width, rects2 + (pos)*RECT_SIZE);
		//return -3;

		out[pos] =  countRectangle(rects1 + pos*RECT_SIZE, rects2 + pos*RECT_SIZE);
	}
}
/*
__global__ void SSIMKernel(unsigned char * data1, unsigned char * data2, double * out, int size, int width){
//
int i = threadIdx.x;
//out[i] =  countSSIM(data1 + size*threadIdx.x, data2 + size*threadIdx.x, size, width);
//countRes(0,0);
}*/



double ** countCUDA(FILE ** streams, FILE * ref, int files_count, PictureData * frame, string type, double ** results) {

	string reference;
	string file1, file2;
	cudaError_t cudaStatus;
	
	const int MAX_BUFFER = 2048000;


	//frame2->data = new char[frame->width*frame->height * 3];
	//cout << frame->frame_count << endl;
	
	//cout << frame->width*frame->height / SKIP_SIZE / SKIP_SIZE*RECT_SIZE << endl;
	//cudaMalloc((void**)&results, frame2->frame_count*sizeof(double));
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		return NULL;
		//goto Error;
	}

	unsigned char ** data=new unsigned char *[files_count];
	unsigned char * dataRef;
	unsigned char ** rects=new unsigned char *[files_count];
	unsigned char ** dataTmp = new unsigned char*[files_count];
	unsigned char * dataTmpTrash = new unsigned char[frame->width*frame->height / 2];
	unsigned char * dataTmpRef = new unsigned char[frame->width*frame->height];
	size_t pitch;

	cudaStatus = cudaMalloc((void **)&dataRef, frame->width*frame->height);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed 1!");
		return NULL;
	}


	//allocated the device memory for source array  
	for (int i = 0; i < files_count; i++) {

		dataTmp[i]= new unsigned char [frame->width*frame->height];
		cudaStatus = cudaMalloc((void **)&data[i], frame->width*frame->height);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed 1!");
			return NULL;
		}
	}
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

	int rec;

	
	for (int i = 0; i < frame->frame_count; i++) {

		rec = fread(dataTmpRef, 1, frame->width*frame->height, ref);
		//fseek(ref, frame->width*frame->height / 2, SEEK_CUR);
		fread(dataTmpTrash, 1, frame->width*frame->height/2, ref); //skip others except Y channel
		if (rec != frame->width*frame->height) {
			printf("error in reading from file 2\n");
			return NULL;
		}
		for (int j = 0; j < files_count; j++) {
			rec = fread(dataTmp[j], 1, frame->width*frame->height, streams[j]);
			if (rec != frame->width*frame->height) {
				printf("error in reading from file 1\n");
				return NULL;
			}
			fread(dataTmpTrash, 1, frame->width*frame->height/2, streams[j]);
			//fseek(streams[j], frame->width*frame->height / 2, SEEK_CUR);//skip others except Y channel

			results[j][i] = countSSIM(dataTmpRef, dataTmp[j], dataRef,data[j], rects[0], rects[1], frame->width*frame->height, frame->width, resultsFrame);
		}
	}

	/*for (int j = 0; j < frame2->frame_count % CHUNK_SIZE; j++) {
	int rec1 = fread(data1[j], 1, frame->width*frame->height * 3, stream);
	int rec2 = fread(data2[j], 1, frame->width*frame->height * 3, stream2);
	if (string(type) == string("SSIM"))  results[frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j] = countSSIM(data1[j], data2[j], frame->width*frame->height, frame->width);
	else results[frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j] = countPSNR(data1[j], data2[j], frame->width*frame->height);
	//cout << frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j << " " << results[frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j] << endl;
	}*/



	//}
	//double * results2=new double[frame2->frame_count];
	/*cudaStatus = cudaMemcpy(results2, results, frame2->frame_count, cudaMemcpyDeviceToHost);*/
	return results;
}



#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>  
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "main.h"   
//#include "..\..\Video_comparsion\Video_comparsion\PSNR.h"   
#include <omp.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <errno.h>
#include <stdio.h>
using namespace std;
//cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);


__device__ float countAvg(unsigned char * data);
__device__ float countVariance(unsigned char * data, float avg);
__device__ float countCovariance(unsigned char * data1, unsigned char * data2, float avg1, float avg2);
__device__ void getRect(unsigned char* data, int start, int width, unsigned char * out);
__device__ float countRectangle(unsigned char * data1, unsigned char * data2);
float countRes(float * tmpRes, int count);
void getLuma(unsigned char *in, unsigned char *out, int size);

float countSSIM(unsigned char * datain1, unsigned char * datain2,unsigned char * dataC1, unsigned char * dataC2, unsigned char * rects1,unsigned char * rects2,int size, int width,float*& results) {
	//unsigned char * data1 = (unsigned char*)datain1;
	//unsigned char * data2 = (unsigned char*)datain2;
	cudaError_t cudaStatus;
	//float * tmpRes = new float[size];
	
	unsigned char * data1 = new unsigned char[size];
	unsigned char * data2 = new unsigned char[size];
/*
	if (data1==0 or data2==0){
		//return -1;
		//cout<<"error in allocation"<<endl;
	}*/
	getLuma(datain1, data1, size);
	getLuma(datain2, data2, size);
	/*unsigned char * rect1 = new unsigned char[RECT_SIZE];
	unsigned char * rect2 = new unsigned char[RECT_SIZE];
	int k = 0;*/

	
	
	//#pragma omp parallel
	//nthreads = omp_get_num_threads();

	
	cudaStatus = cudaMemcpy(dataC1,data1,size, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess ) {
		printf("%s\n", cudaGetErrorString(cudaStatus));
		return -1;
	}
	cudaStatus = cudaMemcpy(dataC2,data2,size, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		printf("%s\n", cudaGetErrorString(cudaStatus));
		return -1;
	}


	//	#pragma omp parallel for schedule(static, 20)
	int rectCount=size/SKIP_SIZE/SKIP_SIZE/THREADS*THREADS;
	int blocks=rectCount/THREADS;
	if (rectCount<size/SKIP_SIZE/SKIP_SIZE){
		blocks=rectCount/THREADS+1;
	}
	countRectangleKernel<<<blocks,THREADS>>>(dataC1,dataC2,rects1,rects2,results,size,width); //FIXME - need to adjust size to count up to THREADS last rectangles!!
	cudaDeviceSynchronize();
	/*for (int i = 0; i < size / width - RECT_SQRT; i += SKIP_SIZE) {

		for (int j = 0; j < width - RECT_SQRT; j += SKIP_SIZE, k++) {
			//for (int i = 0; i < size-(RECT_SQRT-1)*width; i+=SKIP_SIZE) {

			//if (tmpRes[k] < 0) cout << "low result: " << i<< ": " << j<< " :" << tmpRes[k] << endl;

		}
	}


	float res = countRes(tmpRes, k);
	delete[] tmpRes;
	delete[] data1;
	delete[] data2;
	delete[] rect1;
	delete[] rect2;*/
	
	float * resultsOut=new float[rectCount];
	cudaStatus = cudaMemcpy(resultsOut,results,rectCount*sizeof(float), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		printf("%s\n", cudaGetErrorString(cudaStatus));
		return -1;
	}
	float output=countRes(resultsOut,rectCount);
	cout<<output<<endl;
	return output;
}



using namespace std;
void getLuma(unsigned char *in, unsigned char *out, int size) {
	//#pragma omp parallel for schedule(static, 100)
	for (int i = 0; i < size; i++) {
		//out[i] = round(0.216*in[3 * i] + 0.7152*in[3 * i + 1] + 0.0722*in[3 * i + 2]); //get Luma from RGB picture
		//out[i] = round((float)1/3*in[3 * i] + (float)1/3*in[3 * i + 1] + (float)1/3*in[3 * i + 2]); //get Luma from RGB picture
		out[i] = round(0.299*in[3 * i] + 0.587*in[3 * i + 1] + 0.114*in[3 * i + 2]); //get Luma from RGB picture
	}
}




//return one rectangle with RECT_SIZE pixels
__device__ void getRect(unsigned char* data, int start, int width, unsigned char * out) {

	for (int i = 0; i<RECT_SQRT; i++) {
		for (int j = 0; j < RECT_SQRT; j++) {
			out[i*RECT_SQRT + j] = data[start + i*width + j];
		}
		//cudaMemcpy(out + i*RECT_SQRT, data + start + i*width, RECT_SQRT, cudaMemcpyDeviceToDevice);
	}
	//return out;
}

//count ssim of one rectangle with RECT_SIZE pixels
__device__ float countRectangle(unsigned char * data1, unsigned char * data2) {

	float avg1 = countAvg(data1);
	float avg2 = countAvg(data2);

	float var1 = countVariance(data1, avg1);
	float var2 = countVariance(data2, avg2);

	float cov = countCovariance(data1, data2, avg1, avg2);


	float ssim = ((2 * avg1*avg2 + C1)*(2 * cov + C2)) / ((avg1*avg1 + avg2*avg2 + C1)*(var1 + var2 + C2));
	return ssim;
}
//count avg value of given rectangle 
__device__ float countAvg(unsigned char * data) {
	float avg = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		avg += data[i];
	}
	avg = avg / (float)RECT_SIZE;
	return avg;
}

//count variance of given rectangle
__device__ float countVariance(unsigned char * data, float avg) {
	float var = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		var += (data[i] - avg)*(data[i] - avg);
	}
	var = var / (float)RECT_SIZE;
	return var;
}

//count covariance of given rectangle
__device__ float countCovariance(unsigned char * data1, unsigned char * data2, float avg1, float avg2) {
	float cov = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		cov += (data1[i] - avg1)*(data2[i] - avg2);
	}
	cov = cov / (float)RECT_SIZE;
	//if (cov < 0) cout << "neg "<<cov << endl;
	return cov;
}




//count average SSIM value from SSIM values per rectangle
float countRes(float * tmpRes, int count) {
	float sum = 0;
	for (int i = 0; i < count; i += 1) {
		sum += tmpRes[i];

	}
	return sum / (float)count;

}
__global__ void countRectangleKernel(unsigned char * data1, unsigned char * data2,unsigned char * rects1,unsigned char * rects2,float * out,int size, int width){
			int i = threadIdx.x;
			int j= blockIdx.x;
int a;
			int pos=j*THREADS+i;
			if (pos<size/SKIP_SIZE/SKIP_SIZE){
				getRect(data1, ((pos)* SKIP_SIZE)%width + (pos*SKIP_SIZE)/width*SKIP_SIZE*width, width, rects1+(pos)*RECT_SIZE);
				getRect(data2, ((pos)* SKIP_SIZE)%width + (pos*SKIP_SIZE)/width*SKIP_SIZE*width, width, rects2+(pos)*RECT_SIZE);
				//return -3;

				out[pos] = countRectangle(rects1+pos*RECT_SIZE, rects2+pos*RECT_SIZE);
			}
	else{
	rects1[0]=1;
}
a=a+1;
}
/*
__global__ void SSIMKernel(unsigned char * data1, unsigned char * data2, float * out, int size, int width){
	//
	int i = threadIdx.x;
	//out[i] =  countSSIM(data1 + size*threadIdx.x, data2 + size*threadIdx.x, size, width);
	//countRes(0,0);
}*/


int compare(const void * a, const void * b)
{
	return (*(float*)a - *(float*)b);
}
PictureData *getVideoInfo(string path) {
	PictureData * data = new PictureData;
	cout << path.c_str() << endl;
	string cmd = "ffprobe -v error -of flat=s=_  -select_streams v:0 -show_entries stream=width,height,r_frame_rate -show_entries format=duration,nb_frames -of default=noprint_wrappers=1:nokey=1 " + path;
	//string cmd="ffprobe -v error -of flat=s=_ -select_streams v:0 -show_entries stream=width,height,nb_frames -of default=noprint_wrappers=1:nokey=1 "+path;
	string cmd2 = "ffprobe - select_streams v - show_streams" + path + " 2> NUL";

#ifdef __linux__
	FILE *stream = popen(cmd.c_str(), "r");
#else 
	FILE *stream = _popen(cmd.c_str(), "r");
#endif
	char buffer[50];
	fgets(buffer, 10, stream);
	data->width = atoi(buffer);
	fgets(buffer, 10, stream);
	data->height = atoi(buffer);
	fgets(buffer, 20, stream);
	string tmp = buffer;
	int pos = tmp.find('/');
	int fps1 = atoi(buffer);
	float fps2 = atoi(tmp.substr(pos + 1).c_str());
	float fps = fps1 / fps2;
	cout << fps << endl;
	fgets(buffer, 20, stream);
	//cout << buffer << endl;


	float len = atof(buffer);

	cout << len*fps << endl;
	data->frame_count = len*fps;
	//else data->frame_count = 3121;//181250; // 7100;//3121;//1359;//7192;
	return data;
}
void startFFmpeg(string path, FILE *& stream) {
#ifdef __linux__
	string cmd = "ffmpeg -i " + path + " -f image2pipe -pix_fmt rgb24 -vcodec rawvideo - 2>/dev/null";
	cout << cmd << endl;
	stream = popen(cmd.c_str(), "r");
#else 
	string cmd = "ffmpeg -i " + path + " -f image2pipe -threads 3  -pix_fmt rgb24 -vcodec rawvideo - 2>NUL";
	//-c:v h264_qsv
	stream = _popen(cmd.c_str(), "rb");
#endif
	cout << cmd.c_str() << endl;

	if (stream == NULL)
    		printf ("Error opening file: %s\n",strerror(errno));

	//return stream;
}



int main(int argc, char ** argv){
	string reference;
	string file1, file2;
	string type;
	cudaError_t cudaStatus;
	if (argc < 6) { // Check the value of argc. If not enough parameters have been passed, inform user and exit.
		cout << argc << endl;
		cout << "Usage is -r <reference file> -in1 <first video to compare> -in2 <second video to compare> [-type]\n"; // Inform the user of how to use the program
																													   //std::cin.get();
		exit(0);
	}
	else { // if we got enough parameters...

		std::cout << argv[0];
		for (int i = 1; i < argc; i++) { /* We will iterate over argv[] to get the parameters stored inside.
										 * Note that we're starting on 1 because we don't need to know the
										 * path of the program, which is stored in argv[0] */
			if (i + 1 != argc) // Check that we haven't finished parsing already
				if (string(argv[i]) == string("-r")) {
					// We know the next argument *should* be the filename:
					reference = argv[i + 1];
					//std::cout << reference << endl;
				}
				else if (string(argv[i]) == string("-in1")) {
					file1 = string(argv[i + 1]);
					cout << file1.c_str() << endl;
				}
				else if (string(argv[i]) == string("-in2")) {
					file2 = string(argv[i + 1]);
				}
				else if (string(argv[i]) == string("-type")) {
					type = string(argv[i + 1]);

				}
				else {
					//cout << "Not enough or invalid arguments, please try again.\n";
					//Sleep(2000);
					//exit(0);
				}
				std::cout << argv[i] << " ";
		}
	}

	const int MAX_BUFFER = 2048000;

	PictureData * frame = getVideoInfo(file1);
	frame->data = new char[frame->width*frame->height * 3];

	PictureData * frame2 = getVideoInfo(file2);
	frame2->data = new char[frame->width*frame->height * 3];
	cout << frame->frame_count << endl;
	float * results; 

	FILE * stream=0;
	startFFmpeg(file1, stream) ;

	FILE * stream2=0;
	startFFmpeg(file2, stream2) ;
	if (stream == NULL)
    		printf ("Error opening file: %s\n",strerror(errno));
	if (stream2 == NULL)
    		printf ("Error opening file2: %s\n",strerror(errno));

	results=new float[frame2->frame_count];	
//cudaMalloc((void**)&results, frame2->frame_count*sizeof(float));
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		return -1;
		//goto Error;
	}

	unsigned char * data1;
	unsigned char * data2;
	unsigned char * rects1;
	unsigned char * rects2;
	size_t pitch;
	//allocated the device memory for source array  
	cudaStatus=cudaMalloc((void **)&data2, frame->width*frame->height);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		return -1;
	}
	cudaStatus=cudaMalloc((void **)&data1, frame->width*frame->height);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		return -1;
	}
	cudaStatus=cudaMallocPitch((void **)&rects1,&pitch, RECT_SIZE,frame->width*frame->height/SKIP_SIZE/SKIP_SIZE*RECT_SIZE);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		return -1;
	}
	cudaStatus=cudaMallocPitch((void **)&rects2,&pitch, RECT_SIZE,frame->width*frame->height/SKIP_SIZE/SKIP_SIZE*RECT_SIZE);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		return -1;
	}
	
	float * resultsFrame;
	cudaStatus=cudaMalloc((void**)&resultsFrame, frame->width*frame->height/SKIP_SIZE/SKIP_SIZE*sizeof(float));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");

		return -1;

	}

	/*for (int j = 0; j < CHUNK_SIZE; j++) {
		cudaMalloc((void**)&data1[j], frame->width*frame->height * 3);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return -1;
		}
		cudaMalloc((void**)&data2[j], frame->width*frame->height * 3);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			return -1;
		}
		
*/

		/*data1[j] = new unsigned char[frame->width*frame->height * 3];
		data2[j] = new unsigned char[frame->width*frame->height * 3];
	}*/
	unsigned char * datatmp1 = new unsigned char[frame->width*frame->height * 3];
	unsigned char * datatmp2 = new unsigned char[frame->width*frame->height * 3];
	for (int i = 0; i < frame2->frame_count ; i++) {

		int rec1 = fread(datatmp1, 1, frame->width*frame->height * 3, stream);
		if (rec1 != frame->width*frame->height*3) {
			printf("error in reading from file 1\n");
			return -1;
		}		


		int rec2 = fread(datatmp2, 1, frame->width*frame->height * 3, stream2);
		if (rec2 != frame->width*frame->height*3) {
			printf("error in reading from file 2\n");
			return -1;
		}


		//float countSSIM(unsigned char * datain1, unsigned char * datain2,unsigned char * dataC1, unsigned char * dataC2, unsigned char ** rects1,unsigned char ** rects2,int size, int width
		results[i]=countSSIM(datatmp1,datatmp2,data1, data2, rects1,rects2, frame->width*frame->height, frame->width,resultsFrame);
		/*
		omp_set_num_threads(CHUNK_SIZE);
#pragma omp parallel for 
		for (int j = 0; j < CHUNK_SIZE; j++) {
			if (string(type) == string("SSIM")) results[j + i*CHUNK_SIZE] = countSSIM(data1[j], data2[j], frame->width*frame->height, frame->width);
			else results[j + i*CHUNK_SIZE] = countPSNR(data1[j], data2[j], frame->width*frame->height);
			//cout << j+i * CHUNK_SIZE << " " << results[j+i*CHUNK_SIZE] << endl;
		}*/
	}

	/*for (int j = 0; j < frame2->frame_count % CHUNK_SIZE; j++) {
		int rec1 = fread(data1[j], 1, frame->width*frame->height * 3, stream);
		int rec2 = fread(data2[j], 1, frame->width*frame->height * 3, stream2);
		if (string(type) == string("SSIM"))  results[frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j] = countSSIM(data1[j], data2[j], frame->width*frame->height, frame->width);
		else results[frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j] = countPSNR(data1[j], data2[j], frame->width*frame->height);
		//cout << frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j << " " << results[frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j] << endl;
	}*/



	//}
	//float * results2=new float[frame2->frame_count];
	/*cudaStatus = cudaMemcpy(results2, results, frame2->frame_count, cudaMemcpyDeviceToHost);*/
	float sum = 0;
	int frames = frame2->frame_count;
	for (int i = 0; i<frame2->frame_count; i++) {
		cout << i << " " << results[i] << endl;
		if (isfinite(results[i]))
			sum += results[i];
		else frames--;
	}

	cout << "AVG: " << sum / frames << endl;
	qsort(results, frames, sizeof(float), compare);
	cout << "Median: " << results[frames / 2] << endl;
}



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


__device__ double countAvg(unsigned char * data);
__device__ double countVariance(unsigned char * data, double avg);
__device__ double countCovariance(unsigned char * data1, unsigned char * data2, double avg1, double avg2);
__device__ void getRect(unsigned char* data, int start, int width, unsigned char * out);
__device__ double countRectangle(unsigned char * data1, unsigned char * data2);
float countRes(float * tmpRes, int count);
void getLuma(unsigned char *in, unsigned char *out, int size);

double countSSIM(unsigned char * datain1, unsigned char * datain2,unsigned char * dataC1, unsigned char * dataC2, unsigned char ** rects1,unsigned char ** rects2,int size, int width) {
	//unsigned char * data1 = (unsigned char*)datain1;
	//unsigned char * data2 = (unsigned char*)datain2;
	cudaError_t cudaStatus;
	//double * tmpRes = new double[size];
	
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

	float * results;
	cudaStatus=cudaMalloc((void**)&results, size/SKIP_SIZE/SKIP_SIZE);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");

		return -1;

	}
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
	countRectangleKernel<<<size/SKIP_SIZE/SKIP_SIZE/THREADS,THREADS>>>(dataC1,dataC2,rects1,rects2,results,size,width); //FIXME - need to adjust size to count up to THREADS last rectangles!!
	/*for (int i = 0; i < size / width - RECT_SQRT; i += SKIP_SIZE) {

		for (int j = 0; j < width - RECT_SQRT; j += SKIP_SIZE, k++) {
			//for (int i = 0; i < size-(RECT_SQRT-1)*width; i+=SKIP_SIZE) {

			//if (tmpRes[k] < 0) cout << "low result: " << i<< ": " << j<< " :" << tmpRes[k] << endl;

		}
	}


	double res = countRes(tmpRes, k);
	delete[] tmpRes;
	delete[] data1;
	delete[] data2;
	delete[] rect1;
	delete[] rect2;*/

	float * resultsOut=new float[size/SKIP_SIZE/SKIP_SIZE];
	cudaStatus = cudaMemcpy(resultsOut,results,size/SKIP_SIZE/SKIP_SIZE, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		printf("%s\n", cudaGetErrorString(cudaStatus));
		return -1;
	}
	double output=countRes(resultsOut, size/SKIP_SIZE/SKIP_SIZE);
	return output;
}



using namespace std;
void getLuma(unsigned char *in, unsigned char *out, int size) {
	//#pragma omp parallel for schedule(static, 100)
	for (int i = 0; i < size; i++) {
		//out[i] = round(0.216*in[3 * i] + 0.7152*in[3 * i + 1] + 0.0722*in[3 * i + 2]); //get Luma from RGB picture
		//out[i] = round((double)1/3*in[3 * i] + (double)1/3*in[3 * i + 1] + (double)1/3*in[3 * i + 2]); //get Luma from RGB picture
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
__device__ double countRectangle(unsigned char * data1, unsigned char * data2) {

	double avg1 = countAvg(data1);
	double avg2 = countAvg(data2);

	double var1 = countVariance(data1, avg1);
	double var2 = countVariance(data2, avg2);

	double cov = countCovariance(data1, data2, avg1, avg2);


	double ssim = ((2 * avg1*avg2 + C1)*(2 * cov + C2)) / ((avg1*avg1 + avg2*avg2 + C1)*(var1 + var2 + C2));
	return ssim;
}
//count avg value of given rectangle 
__device__ double countAvg(unsigned char * data) {
	double avg = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		avg += data[i];
	}
	avg = avg / (double)RECT_SIZE;
	return avg;
}

//count variance of given rectangle
__device__ double countVariance(unsigned char * data, double avg) {
	double var = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		var += (data[i] - avg)*(data[i] - avg);
	}
	var = var / (double)RECT_SIZE;
	return var;
}

//count covariance of given rectangle
__device__ double countCovariance(unsigned char * data1, unsigned char * data2, double avg1, double avg2) {
	double cov = 0;
	for (int i = 0; i < RECT_SIZE; i++) {
		cov += (data1[i] - avg1)*(data2[i] - avg2);
	}
	cov = cov / (double)RECT_SIZE;
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
__global__ void countRectangleKernel(unsigned char * data1, unsigned char * data2,unsigned char **rects1,unsigned char ** rects2,float * out,int size, int width){
			int i = threadIdx.x;
			int j= blockIdx.x;
			getRect(data1, (j*THREADS+i)* SKIP_SIZE * SKIP_SIZE, width, rects1[i]);
			getRect(data2, (j*THREADS+i)* SKIP_SIZE* SKIP_SIZE, width, rects2[i]);
			//return -3;

			out[i] = countRectangle(rects1[i], rects2[i]);
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
	return (*(double*)a - *(double*)b);
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
	double fps2 = atoi(tmp.substr(pos + 1).c_str());
	double fps = fps1 / fps2;
	cout << fps << endl;
	fgets(buffer, 20, stream);
	//cout << buffer << endl;


	double len = atof(buffer);

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
	unsigned char ** rects1;
	unsigned char ** rects2;
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
	cudaStatus=cudaMallocPitch((void **)&rects1,&pitch, RECT_SIZE,frame->width*frame->height/SKIP_SIZE/SKIP_SIZE);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		return -1;
	}
	cudaStatus=cudaMallocPitch((void **)&rects2,&pitch, RECT_SIZE,frame->width*frame->height/SKIP_SIZE/SKIP_SIZE);
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


		//double countSSIM(unsigned char * datain1, unsigned char * datain2,unsigned char * dataC1, unsigned char * dataC2, unsigned char ** rects1,unsigned char ** rects2,int size, int width
		results[i]=countSSIM(datatmp1,datatmp2,data1, data2, rects1,rects2, frame->width*frame->height, frame->width);
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
	double sum = 0;
	int frames = frame2->frame_count;
	for (int i = 0; i<frame2->frame_count; i++) {
		cout << i << " " << results[i] << endl;
		if (isfinite(results[i]))
			sum += results[i];
		else frames--;
	}

	cout << "AVG: " << sum / frames << endl;
	qsort(results, frames, sizeof(double), compare);
	cout << "Median: " << results[frames / 2] << endl;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    //addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>  
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "SSIM.h"   
#include "main.h"   
#include "psnr.h"   
#include "stvssim.h"
#include "ssim.cuh"
#include <omp.h>

using namespace std;
int compare(const void * a, const void * b)
{
	return (*(double*)a - *(double*)b);
}
PictureData *getVideoInfo(string path) {
	PictureData * data = new PictureData;
	cout << path.c_str() << endl;
	string cmd = "ffprobe -v error -of flat=s=_  -select_streams v:0 -show_entries stream=width,height,r_frame_rate -show_entries format=duration,nb_frames -of default=noprint_wrappers=1:nokey=1 " + path;
	cout << cmd.c_str() << endl;
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
	data->size = data->width*data->height;
	//else data->frame_count = 3121;//181250; // 7100;//3121;//1359;//7192;
	return data;
}
FILE * startFFmpeg(string path) {
	FILE *stream;
#ifdef __linux__
	string cmd = "ffmpeg -i " + path + " -f image2pipe -pix_fmt yuv420p -vcodec rawvideo - 2>/dev/null";
	cout << cmd << endl;
	stream = popen(cmd.c_str(), "rb");
#else 
	string cmd = "ffmpeg -i " + path + " -f image2pipe -threads 3  -pix_fmt yuv420p -vcodec rawvideo - 2>NUL";
	//-c:v h264_qsv
	stream = _popen(cmd.c_str(), "rb");
#endif
	cout << cmd.c_str() << endl;


	return stream;
}

void shiftData(unsigned char ** data, int size) {
	for (int i = 0; i < FRAME_CNT / 2 + 1; i++) {
		memcpy(data[i], data[i + FRAME_CNT / 2], size);
	}

}

double ** countMetric(FILE ** streams, FILE * ref, int files_count, PictureData * frame, string type, double ** results) {

	//double ** results = new double*[files_count];
	unsigned char *** data = new unsigned char **[files_count];
	unsigned char ** dataRef = new unsigned char *[CHUNK_SIZE];
	unsigned char * dataTrash = new unsigned char[frame->width*frame->height / 2];
	int frames = frame->frame_count;
	int rec;
	for (int k = 0; k < files_count; k++) {
		data[k] = new unsigned char *[CHUNK_SIZE];
		//results[k] = new double [frame->frame_count];
		for (int j = 0; j < CHUNK_SIZE; j++) {
			data[k][j] = new unsigned char[frame->width*frame->height];

		}
	}
	for (int j = 0; j < CHUNK_SIZE; j++) {
		dataRef[j] = new unsigned char[frame->width*frame->height];
	}

	for (int i = 0; i < frame->frame_count / CHUNK_SIZE; i++) {
		for (int j = 0; j < CHUNK_SIZE; j++) {
			rec = fread(dataRef[j], 1, frame->width*frame->height, ref);
			if (rec != frame->width*frame->height) {
				cout << "error" << endl;
				return NULL;
			}
			rec = fread(dataTrash, 1, frame->width*frame->height / 2, ref); //when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate
			if (rec != frame->width*frame->height / 2) {
				cout << "error" << endl;
				return NULL;
			}
			for (int k = 0; k < files_count; k++) {
				int rec = fread(data[k][j], 1, frame->width*frame->height, streams[k]);
				if (rec != frame->width*frame->height) {
					cout << "error" << endl;
					return NULL;
				}
				rec = fread(dataTrash, 1, frame->width*frame->height / 2, streams[k]); //when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate
				if (rec != frame->width*frame->height / 2) {
					cout << "error" << endl;
					return NULL;
				}
			}


		}
		omp_set_num_threads(CHUNK_SIZE);
		for (int k = 0; k < files_count; k++) {
#pragma omp parallel for 
			for (int j = 0; j < CHUNK_SIZE; j++) {

				if (string(type) == string("SSIM")) results[k][j + i*CHUNK_SIZE] = countSSIM(dataRef[j], data[k][j], frame->width*frame->height, frame->width);

				else results[k][j + i*CHUNK_SIZE] = countPSNR(dataRef[j], data[k][j], frame->width*frame->height);
				//cout << j+i * CHUNK_SIZE << " " << results[k][j+i*CHUNK_SIZE] << endl;
			}
		}
	}
	for (int j = 0; j < frame->frame_count % CHUNK_SIZE; j++) {
		rec = fread(dataRef[j], 1, frame->width*frame->height, ref);
		fread(dataTrash, 1, frame->width*frame->height/2, ref);
		
		//fseek(ref, frame->width*frame->height / 2, SEEK_CUR); //skip others except Y channel
		if (rec != frame->width*frame->height/2) {
			cout << "error" << endl;
			return NULL;
		}
	}
	for (int k = 0; k < files_count; k++) {
		for (int j = 0; j < frame->frame_count % CHUNK_SIZE; j++) {
			rec = fread(data[k][j], 1, frame->width*frame->height, streams[k]);
			fread(dataTrash, 1, frame->width*frame->height/2, streams[k]);
			
			//fseek(streams[k], frame->width*frame->height / 2, SEEK_CUR); //skip others except Y channel
			if (rec != frame->width*frame->height) {
				cout << "error" << endl;
				return NULL;
			}

			if (string(type) == string("SSIM"))  results[k][frame->frame_count - frame->frame_count % CHUNK_SIZE + j] = countSSIM(dataRef[j], data[k][j], frame->width*frame->height, frame->width);
			else results[k][frame->frame_count - frame->frame_count % CHUNK_SIZE + j] = countPSNR(dataRef[j], data[k][j], frame->width*frame->height);
			//cout << frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j << " " << results[frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j] << endl;
		}
	}

	//delete frame->data;
	delete frame;
	return results;
}


int main(int argc, char ** argv) {

	string reference;
	string file1, file2;
	string type;
	int gpu = 0;
	string  * files = new string[MAX_FILES];
	int files_count = 0;
	if (argc < 6) { // Check the value of argc. If not enough parameters have been passed, inform user and exit.
		cout << argc << endl;
		cout << "Usage is -r <reference file> -in <first video to compare> -in <second video to compare> [-type]\n"; // Inform the user of how to use the program
																													 //std::cin.get();
		exit(0);
	}
	else { // if we got enough parameters...

		//std::cout << argv[0];
		for (int i = 1; i < argc; i++) { /* We will iterate over argv[] to get the parameters stored inside.
										 * Note that we're starting on 1 because we don't need to know the
										 * path of the program, which is stored in argv[0] */
			if (i + 1 != argc) { // Check that we haven't finished parsing already
				if (string(argv[i]) == string("-r")) {
					// We know the next argument *should* be the filename:
					reference = argv[i + 1];
					//std::cout << reference << endl;
				}
				else if (string(argv[i]) == string("-in")) {
					files[files_count] = string(argv[i + 1]);
					//cout << files[files_count].c_str() << endl;
					files_count++;
				}
				/*else if (string(argv[i]) == string("-in2")) {
				file2 = string(argv[i + 1]);
				}*/
				else if (string(argv[i]) == string("-type")) {
					type = string(argv[i + 1]);

				}
				
				else {
					//cout << "Not enough or invalid arguments, please try again.\n";
					//Sleep(2000);
					//exit(0);
				}
				//std::cout << argv[i] << " ";
			}
			else if (i != argc){
				if (string(argv[i]) == string("CUDA")) { //we will use CUDA computation
					gpu = 1;
				}
			}
		}
	}

	const int MAX_BUFFER = 2048000;
	PictureData * frame;
	FILE ** streams;
	FILE * ref;
	streams = new FILE *[files_count];
	double ** results = new double *[files_count];
	frame = getVideoInfo(reference);
	ref = startFFmpeg(reference);
	for (int i = 0; i < files_count; i++) {
		frame = getVideoInfo(files[i]);
		streams[i] = startFFmpeg(files[i]);
		results[i] = new double[frame->frame_count];
	}


	int rec;
	double *sum = new double[files_count];
	int * frames = new int[files_count];
	for (int i = 0; i < files_count; i++) {
		frames[i] = frame->frame_count;
		sum[i] = 0;
	}
	if (gpu == 1) {
		countCUDA(streams, ref, files_count, frame, type, results);
	}

	else if (string(type) == string("STVSSIM")) {

		countMetricSTVSSIM(streams, ref, files_count, frame, type, results);
		//results= new double *[files_count];
		unsigned char ** ref_data = new unsigned char *[FRAME_CNT];
		unsigned char *** data = new unsigned char **[FRAME_CNT];
		for (int k = 0; k < files_count; k++) {
			//results[k] = new double[frame->frame_count];
			data[k] = new unsigned char *[FRAME_CNT];
			for (int j = 0; j < FRAME_CNT; j++) {
				data[k][j] = new unsigned char[frame->size];
			}
		}
		for (int j = 0; j < FRAME_CNT; j++) {
			ref_data[j] = new unsigned char[frame->size];
		}
		unsigned char * tmp = new  unsigned char[frame->size * 3];

		for (int i = FRAME_CNT / 2; i < FRAME_CNT; i++) {
			for (int j = 0; j < files_count; j++) {
				rec = fread(tmp, 1, frame->width*frame->height * 3, streams[j]);
				if (rec != frame->width*frame->height * 3) {
					cout << "error" << endl;
					return -1;
				}
				//getLuma(tmp, data[j][i], frame->size);
			}

			rec = fread(tmp, 1, frame->width*frame->height * 3, ref);
			if (rec != frame->width*frame->height * 3) {
				cout << "error2" << endl;
				return -1;
			}
			//getLuma(tmp, ref_data[i], frame->size);
		}

		int i = FRAME_SKIP;
		int j = 0;

		for (; i < frame->frame_count - FRAME_SKIP; i += FRAME_SKIP, j++) {
			shiftData(ref_data, frame->size);
			for (int k = 0; k < files_count; k++) {
				shiftData(data[k], frame->size);
			}
			for (int k = 0; k < FRAME_SKIP; k++) {
				for (int l = 0; l < files_count; l++) {
					rec = fread(tmp, 1, frame->width*frame->height * 3, streams[l]);
					if (rec != frame->width*frame->height * 3) {
						cout << "error" << endl;
						return -1;
					}
					//getLuma(tmp, data[l][k], frame->size);
				}
			}
			for (int l = 0; l < files_count; l++) {
				results[l][j] = countSTVSSIM(ref_data, data[l], frame->width*frame->height, frame->width);
				cout << j << ": " << results[l][j] << endl;
			}
			//cout << results[j] << endl;

		}
		for (int i = 0; i < files_count; i++) {
			frames[i] = j;
		}
		
	}
	else {
		countMetric(streams, ref, files_count, frame, type, results); //SSIM, PSNR
		
	}

	for (int j = 0; j < files_count; j++) {
		cout << "input file number: " << j << endl;
		for (int i = 0; i < frame->frame_count; i++) {
			cout << i << " " << results[j][i] << endl;
			if (isfinite(results[j][i]))
				sum[j] += results[j][i];
			else frames[j]--;
		}

	}
	for (int i = 0; i < files_count; i++) {
		cout << "AVG: " << sum[i] / frames[i] << endl;
		qsort(results[i], frames[i], sizeof(double), compare);
		cout << "Median: " << results[i][frames[i] / 2] << endl;
	}
}
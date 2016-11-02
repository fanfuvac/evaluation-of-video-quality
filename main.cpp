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
	string cmd = "ffmpeg -i " + path + " -f image2pipe -pix_fmt rgb24 -vcodec rawvideo - 2>/dev/null";
	cout << cmd << endl;
	stream = popen(cmd.c_str(), "rb");
#else 
	string cmd = "ffmpeg -i " + path + " -f image2pipe -threads 3  -pix_fmt rgb24 -vcodec rawvideo - 2>NUL";
	//-c:v h264_qsv
	stream = _popen(cmd.c_str(), "rb");
#endif
	cout << cmd.c_str() << endl;


	return stream;
}

void shiftData(unsigned char ** data,int size) {
	for (int i = 0; i < FRAME_CNT / 2 + 1; i++) {
		memcpy(data[i], data[i + FRAME_CNT / 2], size);
	}

}


int main(int argc, char ** argv) {
	/*#ifdef LINUX
	string file1 = "/cygdrive/d/Dokumenty/Owncloud/FIT/DP/data/W1.mkv";
	string file2 = "/cygdrive/d/Dokumenty/Owncloud/FIT/DP/data/W6.mp4";
	#else
	string file1 = "D:\\Dokumenty\\Owncloud\\FIT\\DP\\data\\W1.mkv";
	string file2 = "D:\\Dokumenty\\Owncloud\\FIT\\DP\\data\\W2.mkv";
	#endif*/
	//int frame_count=7192;

	/*string file1;
	string file2;*/
	string reference;
	string file1, file2;
	string type;
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
	double * results = new double[frame2->frame_count];
	FILE * stream = startFFmpeg(file1);
	FILE * stream2 = startFFmpeg(file2);

	
	double sum = 0;
	int frames = frame2->frame_count;

	if (string(type) == string("STVSSIM")) {

		unsigned char ** data1 = new unsigned char *[FRAME_CNT];
		unsigned char ** data2 = new unsigned char *[FRAME_CNT];
		for (int j = 0; j < FRAME_CNT; j++) {
			data1[j] = new unsigned char[frame->size];
			data2[j] = new unsigned char[frame->size];
		}
		unsigned char * tmp =new  unsigned char[frame->size*3];
		
		for (int i = FRAME_CNT/2; i < FRAME_CNT; i++) {
			int rec1 = fread(tmp, 1, frame->width*frame->height * 3, stream);
			if (rec1 != frame->width*frame->height * 3) {
				cout << "error" << endl;
				return -1;
			}
			getLuma(tmp, data1[i], frame->size);
			int rec2 = fread(tmp, 1, frame->width*frame->height * 3, stream2);
			if (rec2 != frame->width*frame->height * 3) {
				cout << "error2" << endl;
				return -1;
			}
			getLuma(tmp, data2[i], frame->size);
		}

		int i = FRAME_SKIP;
		int j = 0;

		for (; i < frame2->frame_count- FRAME_SKIP; i+=FRAME_SKIP,j++) {
			shiftData(data2, frame->size);
			shiftData(data1, frame->size);
			for (int k = 0; k < FRAME_SKIP; k++) {
				int rec1 = fread(tmp, 1, frame->width*frame->height * 3, stream);
				if (rec1 != frame->width*frame->height * 3) {
					cout << "error" << endl;
					return -1;
				}
				getLuma(tmp, data1[k], frame->size);

				int rec2 = fread(tmp, 1, frame->width*frame->height * 3, stream2);
				if (rec2 != frame->width*frame->height * 3) {
					cout << "error2" << endl;
					return -1;
				}
				getLuma(tmp, data2[k], frame->size);
			}
			results[j]= countSTVSSIM(data1, data2, frame->width*frame->height, frame->width);
			//cout << results[j] << endl;

		}
		frames = j;
		for (int i = 0; i < j; i++) {
			cout << i << " " << results[i] << endl;
				sum += results[i];
		}
	}
	else {
		unsigned char ** data1 = new unsigned char *[CHUNK_SIZE];
		unsigned char ** data2 = new unsigned char *[CHUNK_SIZE];
		for (int j = 0; j < CHUNK_SIZE; j++) {
			data1[j] = new unsigned char[frame->width*frame->height * 3];
			data2[j] = new unsigned char[frame->width*frame->height * 3];
		}

		for (int i = 0; i < frame2->frame_count / CHUNK_SIZE; i++) {
			for (int j = 0; j < CHUNK_SIZE; j++) {
				int rec1 = fread(data1[j], 1, frame->width*frame->height * 3, stream);
				if (rec1 != frame->width*frame->height * 3) {
					cout << "error" << endl;
					return -1;
				}

				int rec2 = fread(data2[j], 1, frame->width*frame->height * 3, stream2);
				if (rec2 != frame->width*frame->height * 3) {
					cout << "error2" << endl;
					return -1;
				}
			}
			omp_set_num_threads(CHUNK_SIZE);
#pragma omp parallel for 
			for (int j = 0; j < CHUNK_SIZE; j++) {
				if (string(type) == string("SSIM")) results[j + i*CHUNK_SIZE] = countSSIM(data1[j], data2[j], frame->width*frame->height, frame->width);
				
				else results[j + i*CHUNK_SIZE] = countPSNR(data1[j], data2[j], frame->width*frame->height);
				//cout << j+i * CHUNK_SIZE << " " << results[j+i*CHUNK_SIZE] << endl;
			}
		}

		for (int j = 0; j < frame2->frame_count % CHUNK_SIZE; j++) {
			int rec1 = fread(data1[j], 1, frame->width*frame->height * 3, stream);
			int rec2 = fread(data2[j], 1, frame->width*frame->height * 3, stream2);
			if (string(type) == string("SSIM"))  results[frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j] = countSSIM(data1[j], data2[j], frame->width*frame->height, frame->width);
			else results[frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j] = countPSNR(data1[j], data2[j], frame->width*frame->height);
			//cout << frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j << " " << results[frame2->frame_count - frame2->frame_count % CHUNK_SIZE + j] << endl;
		}



		//}

		
		
		for (int i = 0; i < frame2->frame_count; i++) {
			cout << i << " " << results[i] << endl;
			if (isfinite(results[i]))
				sum += results[i];
			else frames--;
		}

		delete frame->data;
		delete frame2->data;
		delete frame;
		delete frame2;
	}
	cout << "AVG: " << sum / frames << endl;
	qsort(results, frames, sizeof(double), compare);
	cout << "Median: " << results[frames / 2] << endl;
}
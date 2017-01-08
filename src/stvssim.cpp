#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>  
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "stvssim.h"
#include "SSIM.h"
#include <stdlib.h>     /* abs */
#include <map>
#include <omp.h>


using namespace std;

///count STVSSIM metric for all files 
double ** countMetricSTVSSIM(FILE ** streams, FILE * ref, int files_count, PictureData * frame, double ** results, int *& frames) {
	int rec;
	unsigned char ** ref_data = new unsigned char *[FRAME_CNT];
	unsigned char *** data = new unsigned char **[FRAME_CNT];
	unsigned char * dataTrash = new unsigned char[frame->width*frame->height / 2];
	for (int k = 0; k < files_count; k++) {
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
			readFromFile(data[j][i], frame->size, streams[j]);
			readFromFile(dataTrash, frame->size / 2, streams[j]);//when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate

		}
		readFromFile(ref_data[i], frame->size, ref);
		readFromFile(dataTrash, frame->size / 2, ref);
	}

	int i = FRAME_SKIP;
	int j = 0;

	for (; i < frame->frame_count - 1*FRAME_SKIP-1; i += FRAME_SKIP, j++) {
		shiftData(ref_data, frame->size);
		for (int k = 0; k < files_count; k++) {
			shiftData(data[k], frame->size);
		}
		for (int k = FRAME_CNT / 2 + 1; k < FRAME_CNT; k++) {
			for (int l = 0; l < files_count; l++) {
				readFromFile(data[l][k], frame->size, streams[l]);
				readFromFile(dataTrash, frame->size / 2, streams[l]);//when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate
			}
			readFromFile(ref_data[k], frame->size, ref);
			readFromFile(dataTrash, frame->size / 2, ref);//when using yuv, first 2/3 of the picture are Lumma, others are UV which we do not evaluate

		}
		double resSSIM, res3D;
		for (int l = 0; l < files_count; l++) {
			res3D = countSTVSSIM(ref_data, data[l], frame->width*frame->height, frame->width);
			resSSIM = countSSIM(ref_data[FRAME_CNT / 2], data[l][FRAME_CNT / 2], frame->size, frame->width);
			results[l][j] = res3D*resSSIM;
			//cout<< l<<"\t"<< res3D << "\t" << resSSIM<< "\t"<<results[l][j] << endl;
			//cout << "3D: " << res3D << " SSIM: " << resSSIM << " Total: " << results[l][j] << endl;
		}
	}
	for (int i = 0; i < files_count; i++) {
		frames[i] = j;
	}
	return results;
}

double countSTVSSIM(unsigned char ** datain1, unsigned char ** datain2, int size, int width) {

	unsigned char ** out = new unsigned char*[CHUNK_SIZE];
	int T = ROOD_SIZE;
	vector vct;
	vct.x = 0;
	vct.y = 0;
	int filter;
	unsigned char ***** filters = new unsigned char ****[CHUNK_SIZE]; //generateFilters();
	int rectCount = (int)((size / width - RECT_SQRT_3D) / SKIP_SIZE + 1)*	(int)((width - RECT_SQRT_3D) / SKIP_SIZE + 1);
	double * tmpRes = new double[rectCount];
	unsigned char **** cube1 = new unsigned char ***[CHUNK_SIZE]; //generateCube();
	unsigned char **** cube2 = new unsigned char ***[CHUNK_SIZE]; //generateCube();
	int k = 0;
	for (int i = 0; i < CHUNK_SIZE; i++) {
		filters[i] = generateFilters();
		cube1[i] = generateCube();
		cube2[i] = generateCube();
		out[i] = new unsigned char[RECT_SIZE];
	}
	omp_set_num_threads(CHUNK_SIZE);
#pragma omp parallel for private(vct, filter,k)
	for (int i = 0; i < size / width - RECT_SQRT_3D; i += SKIP_SIZE) {
		int thr = omp_get_thread_num();

		for (int j = 0; j < width - RECT_SQRT_3D; j += SKIP_SIZE) {
			k = (j / SKIP_SIZE) % (width - RECT_SQRT_3D) + (i / SKIP_SIZE) * ((width - RECT_SQRT_3D)/SKIP_SIZE+1); //actual position in loop
			getRect(datain1[FRAME_CNT / 2], i*width + j, width, out[thr]);
			vct = countARPS(out[thr], datain1[FRAME_CNT / 2 - 1], j, i, width, size / width, T);

			if ((vct.x > vct.y * 2 && vct.x*-1 < 2 * vct.y) || (vct.x < vct.y * 2 && vct.x*-1 > 2 * vct.y)) { //y=0
				filter = 0;
			}
			else if ((vct.y > vct.x * 2 && vct.y > -2 * vct.x) || (vct.y < vct.x * 2 && vct.y < -2 * vct.x)) { //x=0
				filter = 2;
			}
			else if ((vct.y > vct.x / 2 && vct.y < 2 * vct.x) || (vct.y < vct.x / 2 && vct.y > 2 * vct.x)) { //y=x
				filter = 3;
			}
			else if ((vct.y > vct.x / -2 && vct.y < -2 * vct.x) || (vct.y <vct.x / -2 && vct.y>-2 * vct.x)) { //y=-x
				filter = 1;
			}
			else if (vct.x == 0 && vct.y == 0) {
				filter = 8;
			}
			else if (vct.x == vct.y * 2) { //exactly between 2 axes
				filter = 4;
			}
			else if (vct.x == vct.y * -2) {
				filter = 5;
			}
			else if (-2 * vct.x == vct.y) {
				filter = 6;
			}
			else if (2 * vct.x == vct.y) {
				filter = 7;
			}
			else {
				cout << "WUT - nonsense vector " << vct.x << " " << vct.y << endl;
			}
			//printf("vct: %d %d %d %d\n",vct.x,vct.y,k,filter);


			//3D-SSIM part
			fillCube(datain1, i*width + j, cube1[thr], width);
			fillCube(datain2, i*width + j, cube2[thr], width);

			//printf("%d %f %f %f %f %d\n",k, res0,res1,res2,res3, filter);*/
			if (filter < 4) {
				tmpRes[k] = countSSIM3D(filters[thr][filter], cube1[thr], cube2[thr]);
			}
			else if (filter < 8) {
				double a, b;
				switch (filter) {
				case 4:
					a = countSSIM3D(filters[thr][0], cube1[thr], cube2[thr]);
					b = countSSIM3D(filters[thr][3], cube1[thr], cube2[thr]);
					break;
				case 5:
					a = countSSIM3D(filters[thr][0], cube1[thr], cube2[thr]);
					b = countSSIM3D(filters[thr][1], cube1[thr], cube2[thr]);
					break;
				case 6:
					a = countSSIM3D(filters[thr][1], cube1[thr], cube2[thr]);
					b = countSSIM3D(filters[thr][2], cube1[thr], cube2[thr]);
					break;
				case 7:
					a = countSSIM3D(filters[thr][2], cube1[thr], cube2[thr]);
					b = countSSIM3D(filters[thr][3], cube1[thr], cube2[thr]);
					break;
				}
				tmpRes[k] = (a + b) / 2;
			}
			else {
				double a, b, c, d;
				a = countSSIM3D(filters[thr][0], cube1[thr], cube2[thr]);
				b = countSSIM3D(filters[thr][1], cube1[thr], cube2[thr]);
				c = countSSIM3D(filters[thr][2], cube1[thr], cube2[thr]);
				d = countSSIM3D(filters[thr][3], cube1[thr], cube2[thr]);
				tmpRes[k] = (a + b + c + d) / 4;
			}
			if (tmpRes[k] >  1) {
				//cout<< tmpRes[k] <<endl;
			}
		}
	}
	double res = countRes(tmpRes, rectCount);
	delete[] tmpRes;

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

	return res;

}

///generate cube filters, vertical, horizontal and 2 inclined are being created
unsigned char **** generateFilters() {
	unsigned char **** out = new unsigned char***[4];
	for (int i = 0; i < 4; i++) {
		out[i] = new unsigned char**[FRAME_CNT];
		for (int j = 0; j < FRAME_CNT; j++) {
			out[i][j] = new unsigned char*[RECT_SQRT_3D];
			for (int k = 0; k < RECT_SQRT_3D; k++) {
				out[i][j][k] = new unsigned char[RECT_SQRT_3D];
			}
		}
	}
	for (int i = 0; i < FRAME_CNT; i++) { //FIXME modify to create 2D array and memcpy
		for (int j = 0; j < RECT_SQRT_3D; j++) {
			for (int k = 0; k < RECT_SQRT_3D; k++) {
				if (j == RECT_SQRT_3D / 2) { //horizontal filter
					out[0][i][j][k] = 1;
				}
				else {
					out[0][i][j][k] = 0;
				}
				if (k == RECT_SQRT_3D / 2) { //vertical filter
					out[1][i][j][k] = 1;
				}
				else {
					out[1][i][j][k] = 0;
				}


				if (j == k) { //x=y filter
					out[2][i][j][k] = 1;
				}
				else {
					out[2][i][j][k] = 0;
				}
				if (k + j == RECT_SQRT_3D) { //x=-y filter
					out[3][i][j][k] = 1;
				}
				else {
					out[3][i][j][k] = 0;
				}
			}
		}
	}
	return out;
}

double countDeltaSqr(unsigned char*** filter, unsigned char*** cube, double mu) {
	double res = 0;
	for (int alpha = 0; alpha < RECT_SQRT_3D; alpha++) {
		for (int beta = 0; beta < RECT_SQRT_3D; beta++) {
			for (int gamma = 0; gamma < FRAME_CNT; gamma++) {
				res += filter[gamma][alpha][beta] * (cube[gamma][alpha][beta] - mu)*(cube[gamma][alpha][beta] - mu);
			}
		}
	}
	return res / (RECT_SQRT_3D*FRAME_CNT);
}

double countDelta(unsigned char*** filter, unsigned char*** cube1, unsigned char*** cube2, double muX, double muY) {
	double res = 0;
	for (int alpha = 0; alpha < RECT_SQRT_3D; alpha++) {
		for (int beta = 0; beta < RECT_SQRT_3D; beta++) {
			for (int gamma = 0; gamma < FRAME_CNT; gamma++) {
				res += filter[gamma][alpha][beta] * (cube1[gamma][alpha][beta] - muX)*(cube2[gamma][alpha][beta] - muX);
			}
		}
	}
	return res / (RECT_SQRT_3D*FRAME_CNT);
}

///Generates 3D array used for SSIM 3D
unsigned char *** generateCube() {
	unsigned char *** out = new unsigned char**[FRAME_CNT];
	for (int i = 0; i < FRAME_CNT; i++) {
		out[i] = new unsigned char*[RECT_SQRT_3D];
		for (int j = 0; j < RECT_SQRT_3D; j++) {
			out[i][j] = new unsigned char[RECT_SQRT_3D];
		}
	}

	return out;
}

///move (FRAME_CNT / 2 + 1) frames from end to begining of the buffer
void shiftData(unsigned char ** data, int size) {
	for (int i = 0; i < FRAME_CNT / 2 + 1; i++) {
		memcpy(data[i], data[i + FRAME_CNT / 2], size);
	}

}

double countSSIM3D(unsigned char *** filter, unsigned char ***  cube1, unsigned char ***  cube2) {
	double muX = countMu(filter, cube1);
	//printf("%f ", muX);
	double muY = countMu(filter, cube2);

	double deltaSqrX = countDeltaSqr(filter, cube1, muX);
	double deltaSqrY = countDeltaSqr(filter, cube2, muY);

	double delta = countDelta(filter, cube1, cube2, muX, muY);

	double ssim3D = ((2 * muX*muY + C1)*(2 * delta + C2)) / ((muX*muX + muY*muY + C1)*(deltaSqrX + deltaSqrY + C2));
	return ssim3D;
}



double countMu(unsigned char*** filter, unsigned char*** cube) {
	double res = 0;
	int res2 = 0;
	for (int alpha = 0; alpha < RECT_SQRT_3D; alpha++) {
		for (int beta = 0; beta < RECT_SQRT_3D; beta++) {
			for (int gamma = 0; gamma < FRAME_CNT; gamma++) {
				res2 += filter[gamma][alpha][beta];
				res += filter[gamma][alpha][beta] * cube[gamma][alpha][beta];
				//cout<<res<<" ";
			}
		}
	}
	return res / (RECT_SQRT_3D*FRAME_CNT);
}

///Fill 3D array with data of surroundings pixels
void fillCube(unsigned char ** datain, int pos, unsigned char *** out, int width) {
	for (int i = 0; i < FRAME_CNT; i++) {
		for (int j = 0; j < RECT_SQRT_3D; j++) {
			memcpy(out[i][j], datain[i] + pos + width*j, RECT_SQRT_3D);
		}
	}
}

//other functions are defined in CUDA version, stvssim.cu


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
		if (x - T - RECT_SQRT  > 0) {
			getRect(framePrev, (x - T) + y*width, width, out);
			res[3] = countSAD(block, out);
		}
		else res[3] = INT_MAX;
		if (x + T + RECT_SQRT < width) {
			getRect(framePrev, (x + T) + y*width, width, out);
			res[1] = countSAD(block, out);
		}
		else res[1] = INT_MAX;
		if (y + T + RECT_SQRT < height) {
			getRect(framePrev, x + (y + T)*width, width, out);
			res[2] = countSAD(block, out);
		}
		else res[2] = INT_MAX;
		if (y - T - RECT_SQRT  > 0) {
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


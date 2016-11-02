#pragma once
struct PictureData {
	int size;
	char* data;
	int width;
	int height;
	int frame_count;
};



#define RECT_SIZE 64 //number of pixels in rectangle
#define RECT_SQRT 8 //size of side of rectangle
#define C1 6.5025 //(K1*L )**2
#define C2 58.5225 //(K2*L )**2
#define SKIP_SIZE 8 //how many pixels should be skipped when going to next window, should be pow of 2((have to be?)
#define CHUNK_SIZE 70 //how many frames will be executed at once
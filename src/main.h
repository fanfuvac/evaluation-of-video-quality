#pragma once
#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>    
struct PictureData {
	int size;
	char* data;
	int width;
	int height;
	int frame_count;
};

int readFromFile(unsigned char *& data, int count, FILE * file);


#define MAX_FILES 10 //max amount of comapred files


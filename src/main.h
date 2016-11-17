struct PictureData{
	int size;
	char* data;
	int width;
	int height;
	int frame_count;
};


#define MAX_FILES 10 //max amount of comapred files

#define THREADS 256 //how many parallel threads will be started for CUDA kernel

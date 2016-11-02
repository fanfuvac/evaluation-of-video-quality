struct PictureData{
	int size;
	char* data;
	int width;
	int height;
};

PictureData *process_file(PictureData* data);
PictureData *read_png_file(char* file_name);
#include <nan.h>
#include <fstream>  
#include <math.h>
#include <string.h>

struct MRCFileHeader {
	int 							nx;
	int								ny;
	int								nz;
	int								mode;
	int								nxstart;
	int								nystart;
	int								nzstart;	
	int								mx;
	int								my;
	int								mz;
	int								cellx;
	int								celly;
	int								cellz;
	int								cella;
	int								cellb;
	int								cellg;
	int								mapc;	
	int								mapr;
	int								maps;
	int								dmin;
	int								dmax;
	int								dmean;
	int								ispg;
	int								nsymbt;
	int								extra1;
	int								extra2;
	int								extra3;
	int								extra4;
	int								extra5;
	int								extra6;
	int								extra7;
	int								extra8;
	int								extra9;
	int								extra10;
	int								extra11;
	int								extra12;
	int								extra13;
	int								extra14;
	int								extra15;
	int								extra16;
	int								extra17;
	int								extra18;
	int								extra19;
	int								extra20;
	int								extra21;
	int								extra22;
	int								extra23;
	int								extra24;
	int								extra25;
	int								originx;
	int								originy;
	int								originz;
	char							map[4];
	char							machst[4];
	int								rms;
	int								nlabl;
	char							label1[80];
	char							label2[80];
	char							label3[80];
	char							label4[80];
	char							label5[80];
	char							label6[80];
	char							label7[80];
	char							label8[80];
	char							label9[80];
	char							label10[80];
};

void parseHeader(std::ifstream& file, MRCFileHeader* header){
		file.read((char*)header, 1024);	
}

void parseMode2Data(std::ifstream& file, MRCFileHeader* header, float* data, int frameid){
	
	int	datalength;
	
	datalength = header->nx * header->ny * 4;
	
	file.seekg (1024 + header->nsymbt + (datalength * frameid), std::ios::beg);
	file.read ((char*)data, datalength);	
}

void parseMode6Data(std::ifstream& file, MRCFileHeader* header, float* data, int frameid){
	
	int		datalength;
	short*	tmparray;
	
	datalength = header->nx * header->ny * 2;
	tmparray = new short[header->nx * header->ny];
	
	file.seekg (1024 + header->nsymbt + (datalength * frameid), std::ios::beg);
	file.read ((char*)tmparray, datalength);
	
	for(int i = 0; i < header->nx * header->ny; i++){
		data[i] = tmparray[i];
	}
	
	delete [] tmparray;
}

void reduceData(MRCFileHeader* header, float* data, char* pixels){
	
	int 	arraylength;
	float	mean;
	float	sd;
	float	normal;
	int 	incr;
	int		counter;
	
	arraylength = header->nx * header->ny;
	mean = 0;
	sd = 0;
	incr = (header->nx/1000) + 1;
	
	counter = 0;
	for(int i = 0; i < arraylength; i+=incr){
		mean += data[i];
		counter++;
	}
	
	mean /= counter;
	
	for(int i = 0; i < arraylength; i+=incr){
		sd += (data[i] - mean) * (data[i] - mean);
	}
	
	sd /= (counter - 1);
	
	sd = sqrt(sd);
	
	for(int i = 0; i < arraylength; i++){
		normal = (12.8 * (data[i] - mean) / sd) + 128;
		if (normal > 255){
			pixels[i] = 255;
		}else if(normal < 0){
			pixels[i] = 0;
		}else{
			pixels[i] = normal;
		}
	}
	
}

NAN_METHOD(getHeader) {
	
	std::ifstream 			file;
	MRCFileHeader*			header;
	std::string 			filename;
	v8::Local<v8::Object>	json;
	
	
	header = new MRCFileHeader();
	json = Nan::New<v8::Object>();
	filename = std::string(*Nan::Utf8String(info[0]));
	
	file.open(filename, std::ios::binary);
	
	if(file.is_open()){
		parseHeader(file, header);
		file.close();
	}
	
	Nan::Set(json, Nan::New("nx").ToLocalChecked(), Nan::New(header->nx));
	Nan::Set(json, Nan::New("ny").ToLocalChecked(), Nan::New(header->ny));
	Nan::Set(json, Nan::New("nz").ToLocalChecked(), Nan::New(header->nz));
	Nan::Set(json, Nan::New("mode").ToLocalChecked(), Nan::New(header->mode));
	Nan::Set(json, Nan::New("nsymbt").ToLocalChecked(), Nan::New(header->nsymbt));
	
	info.GetReturnValue().Set(json);
	
	delete header;
	
}

NAN_METHOD(getData) {
	
	std::ifstream 	file;
	MRCFileHeader*	header;
	float*			data;
	
	header = new MRCFileHeader();
	
	file.open("testmic.mrc", std::ios::binary);
	
	if(file.is_open()){
		parseHeader(file, header);
		data = new float[header->nx * header->ny];
		if(header->mode == 2){
			parseMode2Data(file, header, data, 0);
		}
		file.close();
		delete [] data;
	}
	
	auto message = Nan::New("header->nx").ToLocalChecked();
	info.GetReturnValue().Set(message);
	
	delete header;
	
}

NAN_METHOD(toPixels){
	
	std::ifstream 			file;
	MRCFileHeader*			header;
	float*					data;
	char*					pixels;
	std::string 			filename;
	v8::Local<v8::Object>	json;
	
	header = new MRCFileHeader();
	json = Nan::New<v8::Object>();
	filename = std::string(*Nan::Utf8String(info[0]));
	
	file.open(filename, std::ios::binary);
	
	if(file.is_open()){
		parseHeader(file, header);
		data = new float[header->nx * header->ny];
		pixels = new char[header->nx * header->ny];
		if(header->mode == 2){
			parseMode2Data(file, header, data, (int)info[1]->IntegerValue());
		} else if(header->mode == 6){
			parseMode6Data(file, header, data, (int)info[1]->IntegerValue());
		}
		file.close();
		reduceData(header, data, pixels);
		Nan::Set(json, Nan::New("pixbuf").ToLocalChecked(), Nan::CopyBuffer(pixels, header->nx * header->ny).ToLocalChecked());
		delete [] data;	
		delete [] pixels;
	}
	Nan::Set(json, Nan::New("nx").ToLocalChecked(), Nan::New(header->nx));
	Nan::Set(json, Nan::New("ny").ToLocalChecked(), Nan::New(header->ny));
	Nan::Set(json, Nan::New("mode").ToLocalChecked(), Nan::New(header->mode));
	info.GetReturnValue().Set(json);
	delete header;
}

NAN_MODULE_INIT(Initialize) {
    NAN_EXPORT(target, getHeader);
    NAN_EXPORT(target, getData);
    NAN_EXPORT(target, toPixels);
}

NODE_MODULE(addon, Initialize);


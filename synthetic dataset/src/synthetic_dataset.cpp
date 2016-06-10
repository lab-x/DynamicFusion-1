
#include <iostream>
#include <cvd/image_io.h>
#include <cvd/image_ref.h>

using namespace std;

int main() {
	int count = 0;
	CVD::Image<u_int16_t>depthmap(CVD::ImageRef(320 * 2,240 * 2));

//	string inputfileNamepng = "../ball.png";
//	std::cout<<"Trying to read file from: " << inputfileNamepng << std::endl;
//	CVD::Image<u_int16_t>DEPTHPNG;
//	CVD::img_load(DEPTHPNG, inputfileNamepng);
	float scale = 1;
	for (int k = 0; k < 50; k++) {
		for (int i = 0; i < 240*2; i++) {
			for (int j = 0; j < 320*2; j++) {
				depthmap[CVD::ImageRef(j,i)] = 0;
			}
		}
//
//     standing triangle
//		for(int i = 50; i < 100; i++)
//		{
//			for(int j = 50; j < 150; j++)
//			{
//				depthmap[CVD::ImageRef(j,i)] = (7000 + 50 * (j - 50)) * scale;
//			}
//		}
//
//		for(int i = 50; i < 100; i++)
//		{
//			for(int j = 150; j < 250; j++)
//			{
//				depthmap[CVD::ImageRef(j,i)] = (7000 + 50 * (250 - j)) * scale;
//			}
//		}
//
//		scale = scale * 1.01;

//		std::cout<< "read: " << std::endl;
//
//		int h_width = 320;
//		int h_height = 240;
//
//		for (int v = 0; v < h_height; v++)
//		{
//			for (int u = 0; u < h_width; u++)
//			{
//				//if (((float)DEPTHPNG[CVD::ImageRef(u,v)])/5000 < 2)
//					depthmap[CVD::ImageRef(u, v)] = ((float)DEPTHPNG[CVD::ImageRef(u,v)]) * scale;
//			}
//		}
//		scale = scale *1.05;


		for (int i = 0; i < 200*2; i++) {
			int j_s = i > 140*2 ? i - 140*2 : 0;
			int j_e = i < 60 ? i : 60;
			for (int j = j_s; j < j_e; j++) {
				depthmap[CVD::ImageRef(j + 20, i + 20)] = (3 * 5000 + j * 70) * scale;
			}
		}
		for (int i = 0; i < 200-60; i++) {
			for (int j = 0; j < 80; j++) {
				depthmap[CVD::ImageRef(j + 80, i + 80)] = (3 * 5000 + 59 * 70) * scale;
			}
		}
		for (int i = 0; i < 200; i++) {
			int j_s = i < 60 ? 60 - i : 0;
			int j_e = i > 140 ? 200 - i : 60;
			for (int j = j_s; j < j_e; j++) {
				depthmap[CVD::ImageRef(j + 160, i + 20)] = (3 * 5000 + (60-j) * 70) * scale;
			}
		}

		scale = scale * 1.01;

		char fileName[200];
		std::sprintf(fileName,"test_%05d.png",count++);
		CVD::img_save(depthmap, fileName);
	}


	return 0;
}

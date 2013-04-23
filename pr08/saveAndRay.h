#include "cyPoint.h"
using namespace cy;
struct RGBType {
	double r;
	double g;
	double b;
};


void savebmp(const char *filename, int w, int h, int dpi, RGBType *data){
	FILE *f;
	int k = w*h;
	int s = 4*k;
	int filesize = 54+s;
	
	double factor = 39.375;
	int m = static_cast<int>(factor);
	
	int ppm = dpi*m;
	
	unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0,0,0, 54,0,0,0};
	unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
	
	bmpfileheader[ 2] = (unsigned char)(filesize);
	bmpfileheader[ 3] = (unsigned char)(filesize>>8);
	bmpfileheader[ 4] = (unsigned char)(filesize>>16);
	bmpfileheader[ 5] = (unsigned char)(filesize>>24);
	
	bmpinfoheader[ 4] = (unsigned char)(w);
	bmpinfoheader[ 5] = (unsigned char)(w>>8);
	bmpinfoheader[ 6] = (unsigned char)(w>>16);
	bmpinfoheader[ 7] = (unsigned char)(w>>24);
	
	bmpinfoheader[ 8] = (unsigned char)(h);
	bmpinfoheader[ 9] = (unsigned char)(h>>8);
	bmpinfoheader[10] = (unsigned char)(h>>16);
	bmpinfoheader[11] = (unsigned char)(h>>24);
	
	bmpinfoheader[21] = (unsigned char)(s);
	bmpinfoheader[22] = (unsigned char)(s>>8);
	bmpinfoheader[23] = (unsigned char)(s>>16);
	bmpinfoheader[24] = (unsigned char)(s>>24);
	
	bmpinfoheader[25] = (unsigned char)(ppm);
	bmpinfoheader[26] = (unsigned char)(ppm>>8);
	bmpinfoheader[27] = (unsigned char)(ppm>>16);
	bmpinfoheader[28] = (unsigned char)(ppm>>24);
	
	bmpinfoheader[29] = (unsigned char)(ppm);
	bmpinfoheader[30] = (unsigned char)(ppm>>8);
	bmpinfoheader[31] = (unsigned char)(ppm>>16);
	bmpinfoheader[32] = (unsigned char)(ppm>>24);
	
	f = fopen(filename, "wb");
	
	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);
	
	for(int i = 0; i < k; i++){
		RGBType rgb = data[i];
		
		double red 	= (data[i].r)*255;
		double green 	= (data[i].g)*255;
		double blue 	= (data[i].b)*255;
		
		unsigned char color[3] = { (int)floor(blue), (int)floor(green), (int)floor(red)};
		
		fwrite(color,1,3,f);
	}
	
	fclose(f);

};
class Color
{
	public:
	double red,green,blue;
	Color(){red=0;blue=0;green=0;}
	Color(Point3f a)
	{
		red = a.x; green = a.y; blue = a.z;
	}
	Color(double r,double g, double b){
		red = r;green=g;blue=b;
	}
	Color operator+(Color a)
	{
		return Color(red+a.red , green+a.green , blue+a.blue);
	}
	Color operator-(Color a)
	{
		return Color(red-a.red , green-a.green , blue-a.blue);
	}
	Color operator*(Color a)
	{
		return Color(red*a.red , green*a.green , blue*a.blue);
	}
	//*
	void printColor()
	{
		printf("(%g , %g , %g)",red,green,blue);
	}

	Color operator*(float c){
	      return Color(red*c, green*c, blue*c );}
	Color operator*(double c){
		return Color(red*c, green*c, blue*c );}
	Color operator/(double c){
		return Color(red/c, green/c, blue/c );}
        Color operator/(float c){
          return Color(red/c, green/c, blue/c );}
// */
	                                                                      
};	
/*************************************** RAY CLASS **********************************************************************/
class Ray
{
	public:
	Point3f origin,direction,color;
	Ray()
	{	origin=Point3f(0,0,0),direction=Point3f(0,0,0),color=Point3f(0,0,0);	}
	Ray(Point3f rs,Point3f rd)
	{	origin=rs,direction=rd,color=Point3f(1,1,1);direction.Normalize();	}
	Ray(Point3f rs,Point3f rd, Point3f rc)
	{	origin=rs,direction=rd,color=rc;direction.Normalize();			}
};

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "Vect.h"
#include "Ray.h"
#include "Camera.h"
#include "Color.h"
#include "Light.h"
#include "Object.h" // obj class - any obj in scene is a subclass of obj class
#include "Sphere.h"
#include "Plane.h"

using namespace std;

//Create a data type for RGB

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

}

int winningObjectIndex(vector<double> object_intersections){
	//returns the index of winning intersections
	int index_of_minimum_value;
	
	//prevent unnecassarry calculations
	if (object_intersections.size() == 0){
		//if there are no intersections
		return -1;	
	}
	
	else if (object_intersections.size() == 1){
		if (object_intersections.at(0) > 0){
			// if that intersections is greater than 0 then its our index of minimum value
			return 0; // 0 is index value, 0th element
		}
		else{
			//otherwise the only intersection is negative
			return -1;
		}
	}
		else{
			// more than one intersections
			// firest find the max value in vector
			
			double max = 0;
			// object_intersections.size() is the no xsections in our vector
			for (int i =0; i < object_intersections.size(); i++){
				if (max < object_intersections.at(i)){
					max = object_intersections.at(i); //so we go thru intersections and the first one we encounter is greater than 0 its max
					
				}		
			}
			
			// then starting from max valur find the minimum positive value
			if (max > 0){
				// we only want positive intersections
				for (int index = 0; index < object_intersections.size(); index++){
					if(object_intersections.at(index) > 0 && object_intersections.at(index) <= max){
						max = object_intersections.at(index); //reset max to new value
						index_of_minimum_value = index;		//smallest positive number to return
					}
				}
				return index_of_minimum_value;	
			}
			else{
				// all intersections were negative
				return -1;	
			}
		}
	
}

// Main function
int thisone;

int main(int argc, char *argv[]){
	
	
	cout<<"Rendering..\n";

	int dpi = 72;
	
	int width = 640;
	int height = 480;
	
	int n = width*height;
	RGBType *pixels = new RGBType[n];
	
	double aspectratio = (double)width/(double)height;

	Vect O (0,0,0);	// corresponds of origin
	
	Vect X (1,0,0);
	Vect Y (0,1,0);
	Vect Z (0,0,1);

	Vect campos (3, 1.5, -4);

	//define the pos dir and rite and down cordinates for camera using vec opeeration
	Vect look_at (0,0,0);
	Vect diff_btw ( campos.getVectX() - look_at.getVectX(), campos.getVectY() - look_at.getVectY(), campos.getVectZ() - look_at.getVectZ());
	Vect camdir = diff_btw.negative().normalize();
	Vect camright = Y.crossProduct(camdir).normalize();
	Vect camdown = camright.crossProduct(camdir);
	
	Camera scene_cam(campos, camdir, camright, camdown);

	//special values for reflectivity and shininess
	Color white_light(1.0,1.0,1.0,0);
	Color white(1.0,1.0,1.0,0);
	Color pretty_green(0.5,1.0,0.5,0.3);
	Color maroon(0.5, 0.25, 0.25, 0);
	Color unknown(0.8, 0.6, 0.35, 0);
	Color gray(0.5,0.5,0.5, 0);
	Color unknw(0.1,0.6,0.8, 0);
	Color black(0.0,0.0,0.0,0);

	//create a light for light source
	Vect light_position (-7,10,-10);
	Light scene_light(light_position, maroon);
	
	// scene objects
	Sphere scene_sphere (O, 1, pretty_green);
	Plane scene_plane (Y, -1.5, gray); // plane facing with y up, dist from origin is -1 down
	
	//we need to put the objects in stack or array so we canindex them
	vector<Object*> scene_objects;
	scene_objects.push_back(dynamic_cast<Object*>(&scene_sphere));
	scene_objects.push_back(dynamic_cast<Object*>(&scene_plane));
	
	double xamnt, yamnt; //camera directions, ray needs to go L/R U/D in img plane of camera
	
	for(int x=0;x<width;x++){
		for(int y=0;y<height;y++){
			thisone = y * width + x; // value at x,y cordinates of individaul pixels
			
			// start with no anti aliasing
			if(width > height) {
				//the image is wider than it is tall
				xamnt = ((x+0.5)/width)*aspectratio - (((width-height)/(double)height)/2);
				yamnt = ((height-y) + 0.5)/height;
			}
			else if (height > width) {
				//the image is taller than it is wide
				xamnt = (x+0.5)/width;
				yamnt = (((height-y) + 0.5)/height)/aspectratio - (((height - width)/(double)width)/2); 
			}
			else{
				//the image is square
				xamnt = (x + 0.5)/width;
				yamnt = ((height-y)+0.5)/height;
			}
			
			// creating rays, origin for all rays is same as origin of camera
			Vect cam_ray_origin = scene_cam.getCameraPosition(); //this is to return camera's origin
			Vect cam_ray_direction = camdir.vectAdd(camright.vectMult(xamnt-0.5).vectAdd(camdown.vectMult(yamnt-0.5))).normalize(); //cam ray dir goes thru xamnt and yamnt
			
			//new insstance of ray
			Ray cam_ray(cam_ray_origin, cam_ray_direction); //goes thru specific x and y
			
			//now send into scene and find the intersectoins
			vector<double> intersections; //vector array to hold xsection values
			
			//since we are insde for loop, at a specific pixel,
			//we need to loop thru each obj in scene and find if ray xsects with any obj in scne
			
			
			
			for (int index=0;index<scene_objects.size();index++){
				//loops thru each obj to find, ask to findIntersection with camray and push that value into xsetion vector
				//we started off by offseting values from the direction camera was pointing
				//creating rays that slightly go in diff directions diff of camera dir
				//if there is xsections put in xsection array
				intersections.push_back(scene_objects.at(index)->findIntersection(cam_ray));
			}
			//we neeed to find which is closest to camera
			int index_of_winning_object = winningObjectIndex(intersections);
			//is a function that sorts the ray of xsections and returns the index of winning objects
			//ie same index in scene objects
			
			// now we have a value for "index_of_winning_object", which is the object that is close 
			// to camera
			//get the color of that individual x and y pixel
			
			if ( index_of_winning_object == -1 ){
				// it missed so its black, bg = black
				pixels[thisone].r = 0;
				pixels[thisone].g = 0;
				pixels[thisone].b = 0;
			}
			else{
				//index corrsponds to object in our scene
				//we wanna get color of that object
				
				Color this_color = scene_objects.at(index_of_winning_object)->getColor();
				pixels[thisone].r = this_color.getColorRed();
				pixels[thisone].g = this_color.getColorGreen();
				pixels[thisone].b = this_color.getColorBlue();
			}
		}
	}
	savebmp("scene.bmp", width, height, dpi, pixels);
	return 0;
}

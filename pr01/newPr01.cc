#include<iostream>
#include<cstdio>
#include<cmath>
#include<vector>
#include<stack>
#include<queue>
#include<map>
#include<sstream>
#include<algorithm>
#include<string>
#include<limits.h>
#include "cyPoint.h"
#include <ctime>

using namespace cy;
using namespace std;
#define print(p) printf("(%f,%f,%f) \n",p.x,p.y,p.z);
#define pb(a) push_back(a)

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
Point3f Pe(0,0,0);	//camera or eye position
Point3f PL(2,2,0);

/*************************************** RAY CLASS **********************************************************************/
class Color
{
	public:
	double red,green,blue;
	Color(){red=0;blue=0;green=0;}
	Color(Point3f a)
	{
		red = a.x; green = a.y; blue = a.z;
	}
	Color(double R,double G, double B){
		red = R;green=G;blue=B;
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

	//Color operator*(float c){
	//	return Color(red*c, green*c, blue*c );}
	Color operator*(double c){
		return Color(red*c, green*c, blue*c );}
	// */

};

class Ray
{
	public:
	Point3f origin,direction; Color color;
	Ray()
	{	origin=Point3f(0,0,0),direction=Point3f(0,0,0),color=Color(Point3f(0,0,0));	}
	Ray(Point3f rs,Point3f rd)
	{	origin=rs,direction=rd,color=Color(Point3f(1,1,1));direction.Normalize();	}
	Ray(Point3f rs,Point3f rd, Color rc)
	{	origin=rs,direction=rd,color=rc;direction.Normalize();			}
};
/*************************************** OBJECT CLASS **********************************************************************/
class Object 
{
	public:
		string objectName;
		Object(){}
		virtual Point3f getIntersectionPoint(Ray ray){return Point3f(0,0,0);}
		virtual bool isEyeOutside(Point3f Pe){return false;}
		virtual Color getColor(){return Color(0,0,0);}
		virtual Color lambertShader(Ray myray){return Color(0,0,0);}
		virtual Color goochShader(Ray ray){return Color(0,0,0);}

		virtual Color phongShader(Ray ray){return Color(0,0,0);}

};
/*************************************** OBJECT CLASS **********************************************************************/
class Plane: public Object
{
	public:
		
		Point3f normalVector;
		Point3f origin;
		Color color;
		Plane(){normalVector = Point3f(0,0,1);origin=Point3f(0,0,0);color=Color(1,0.5,0.2);objectName="Plane";}
		Plane(Point3f r, Point3f orig){	normalVector = r; origin=orig;	objectName="Plane";normalVector.Normalize();	}
		Plane(Point3f r, Point3f orig,Color c){	normalVector = r; origin=orig; color = c; objectName="Plane";	normalVector.Normalize();}
		virtual bool isEyeOutside(Point3f Pe)
		{
			double c=(Pe-origin)%(normalVector);
			
			if(c>0)		//Eye is outside the sphere 
				return true;
			else 	return false;
		}
		virtual Color getColor(){	return color;}
		virtual Point3f getIntersectionPoint(Ray ray)
		{
			Point3f ret(-1,-1,-1);
			//cout<<"assuming that eye is outside the sphere \n";
			double t= (normalVector%(origin-ray.origin))/(normalVector%ray.direction);

			ret = ray.origin + t*ray.direction;
			return ret;
		}	
		virtual Color lambertShader(Ray myray)
		{
			Color spColor ;
			Point3f Ph,Pc, nh,nlh;
			//print(interPoint);cout<<endl;
			double cosTheta=0,c;
			double alpha,s,s0,delta;
			Point3f v,r;

			//for(int i=0; i < allObjects.size(); i++)
			//{
			spColor = color;
			Ph=getIntersectionPoint(myray);

			Color ambient_color(0,0,0),diffused_color,finalColor,specular_color(1,1,1);
			// ========================== Lambert and Gooch =========================================	
			Pc=origin;
			diffused_color = color;
			//nh=(Ph - Pc);nh.Normalize();
			nh = normalVector;

			nlh=(PL - Ph);nlh.Normalize();
			cosTheta = nh % nlh;

			// try different values with c and cosTheta
			c = cosTheta;
			//c = (cosTheta+1.0)/2.0;
			//cout<< "cos: "<<c<<endl;
			if(c<0)	c=0.0;
			//s=pow(s,alpha);
			//cout<< "c: "<<c<<endl;
			// Gooch shading
			finalColor = diffused_color*c + ambient_color*(1-c);
			//finalColor = finalColor * color;
			//finalColor.printColor();cout<<endl;
			//========================================== L & G =======================================

			v=Pe-Ph;v.Normalize();
			r= 2*(v%nh)*nh-v ;r.Normalize();
			s = nlh%r,s0=s,delta=1;
			//cout<<"s: "<<s<<endl;
			alpha=4;
			//s=pow(s,alpha);
			//cout<<"v.n and r.n: "<<v%nh<<" "<<r%nh<<endl;
			//if (abs(v%nh - r%nh) >= 0.001 )	s=1;
			if (s<0.0)	s=0;
			else if(s>1.0)	s=1;
			//cout<<"s: "<<s<<endl;
			s=pow(s,alpha);
			//cout<<"s: "<<s<<endl;
			//if (s<0)	s=0;
			finalColor = finalColor*(1-s) + specular_color*(s);
			return finalColor;

		}
		virtual Color goochShader(Ray myray)
		{
			Color spColor ;
			Point3f Ph,Pc, nh,nlh;
			double cosTheta=0,c;
			double alpha,s,s0,delta;
			Point3f v,r;
			Color ambient_color(0,0,0),diffused_color,finalColor,specular_color(1,1,1);
			//for(int i=0; i < allObjects.size(); i++)
			//{
			spColor = color;
			Ph=getIntersectionPoint(myray);

			// ========================== Lambert and Gooch =========================================	
			diffused_color = color;
			Pc=origin;
			nh=normalVector;
			nlh=(PL - Ph);nlh.Normalize();
			cosTheta = nh % nlh;

			// try different values with c and cosTheta
			c = cosTheta;
			c = (cosTheta+1.0)/2.0;
			if(c<0)	c=0.0;
			finalColor = diffused_color*c + ambient_color*(1-c);
			return finalColor;
		}

		virtual Color phongShader(Ray myray)
		{
			Color spColor ;
			Point3f Ph,Pc, nh,nlh;
			double cosTheta=0,c;
			double alpha,s,s0,delta;
			Point3f v,r;
			Color ambient_color(0,0,0),diffused_color,finalColor,specular_color(1,1,1);

			Ph=getIntersectionPoint(myray);
			nh=normalVector;
			nlh=(PL - Ph);nlh.Normalize();
			finalColor = lambertShader(myray);

			v=Pe-Ph;v.Normalize();
			r= 2*(v%nh)*nh-v ;r.Normalize();
			s = nlh%r,s0=s,delta=1;
			alpha=4;
			if (s<0.0)	s=0;
			else if(s>1.0)	s=1;
			s=pow(s,alpha);
			finalColor = finalColor*(1-s) + specular_color*(s);
			return finalColor;
		}


};

/************************************** SPHERE CLASS ***********************************************************************/
class Sphere : public Object
{
	public:
		int id;
	double radius;
	Point3f center; Color color;
	Sphere()
	{
		radius= 0.0;
		center=Point3f(0,0,0);
		color=Color(0,0,0);
		objectName = "Sphere";
	}
	Sphere(Point3f c, double rad,Color col,int i)
	{
		radius= rad;
		id=i;
		center=c;
		color=col;
		objectName = "Sphere";
	}
	Sphere(Point3f c, double rad, int i)
	{
		id=i;
		radius= rad;
		center=c;
		color=Color(1,1,1);
		objectName = "Sphere";
	}
	
	Point3f getNormal(Point3f point)
	{
		Point3f n=point-center;
		n.Normalize();
		return n;
	}
	virtual bool isEyeOutside(Point3f Pe)
	{
		double c=(center-Pe)%(center-Pe)-(radius*radius);
		
		if(c>0)		//Eye is outside the sphere 
			return true;
		else 	return false;
	}
	virtual Color getColor(){	return color;}
	virtual Point3f getIntersectionPoint(Ray ray)
	{
		Point3f ret;
		//cout<<"assuming that eye is outside the sphere \n";
		
		double c=(center-ray.origin)%(center-ray.origin)-(radius*radius),b=(center-ray.origin)%ray.direction;
		double delta=b*b-c,t1,t2,tmin;
		if(delta<0)	{
			//cout<<"delta<0\n";
			return Point3f(-1,-1,-1);
		}
		//*
		if(b<0)	{
			//cout<<"b<0\n";
			return Point3f(-1,-1,-1);
		}
		// */
		t1=b+sqrt(delta);t2=b-sqrt(delta);
		//cout<<t1<<" , "<<t2<<endl;
		tmin=min(t1,t2);
		//cout<<"min t is : "<<tmin<<endl;
		ret=Point3f((ray.origin+tmin*ray.direction));
		
		return ret;
	}
	virtual Color lambertShader(Ray myray)
	{
		Color spColor ;
		Point3f Ph,Pc, nh,nlh;
		double cosTheta=0,c;
		double alpha,s,s0,delta;
		Point3f v,r;
		Color ambient_color(0,0,0),diffused_color,finalColor,specular_color(1,1,1);
		//for(int i=0; i < allObjects.size(); i++)
		//{
		spColor = color;
		Ph=getIntersectionPoint(myray);

		// ========================== Lambert and Gooch =========================================	
		diffused_color = color;
		Pc=center;
		nh=(Ph - Pc);nh.Normalize();

		nlh=(PL - Ph);nlh.Normalize();
		cosTheta = nh % nlh;

		// try different values with c and cosTheta
		c = cosTheta;
		//c = (cosTheta+1.0)/2.0;
		if(c<0)	c=0.0;
		finalColor = diffused_color*c + ambient_color*(1-c);
		return finalColor;

	}
	virtual Color goochShader(Ray myray)
	{
		Color spColor ;
		Point3f Ph,Pc, nh,nlh;
		double cosTheta=0,c;
		double alpha,s,s0,delta;
		Point3f v,r;
		Color ambient_color(0,0,0),diffused_color,finalColor,specular_color(1,1,1);
		//for(int i=0; i < allObjects.size(); i++)
		//{
		spColor = color;
		Ph=getIntersectionPoint(myray);

		// ========================== Lambert and Gooch =========================================	
		diffused_color = color;
		Pc=center;
		nh=(Ph - Pc);nh.Normalize();

		nlh=(PL - Ph);nlh.Normalize();
		cosTheta = nh % nlh;

		// try different values with c and cosTheta
		c = cosTheta;
		c = (cosTheta+1.0)/2.0;
		if(c<0)	c=0.0;
		finalColor = diffused_color*c + ambient_color*(1-c);
		return finalColor;
	}

	virtual Color phongShader(Ray myray)
	{
		//========================================== L & G =======================================
		Color spColor ;
		Point3f Ph,Pc, nh,nlh;
		double cosTheta=0,c;
		double alpha,s,s0,delta;
		Point3f v,r;
		Color ambient_color(0,0,0),diffused_color,finalColor,specular_color(1,1,1);

		finalColor = lambertShader(myray);
		spColor = color;
		Ph=getIntersectionPoint(myray);

		// ========================== Lambert and Gooch =========================================	
		diffused_color = color;
		Pc=center;
		nh=(Ph - Pc);nh.Normalize();

		nlh=(PL - Ph);nlh.Normalize();
		
		v=Pe-Ph;v.Normalize();
		r= 2*(v%nh)*nh-v ;r.Normalize();
		s = nlh%r,s0=s,delta=1;
		alpha=4;
		if (s<0.0)	s=0;
		else if(s>1.0)	s=1;
		s=pow(s,alpha);
		finalColor = finalColor*(1-s) + specular_color*(s);
		return finalColor;
	}

};
int findWinningPointIndex(vector<Point3f > myinter ,Point3f Pe)
{
	float len=0,minLen=INT_MAX,index=-1;
	for (int i = 0; i < myinter.size(); i++)
	{
		if(myinter[i].x == -1 && myinter[i].y == -1 && myinter[i].z == -1)	continue;
		len=(myinter[i]-Pe).Length();
		if(len < minLen)
		{
			minLen=len;
			index=i;
		}
	}
	
	return index;
}

/********************************************* MAIN ****************************************************************/
int main (int argc, char const* argv[])
{
	int antiAliasing=0;
	
	
	int Xmax = 500;
	int Ymax = 500;
	int n=Ymax*Xmax,dpi=72;
	
	
	RGBType *pixels= new RGBType[n];
	int index=0,M=2,N=2; 
	int Sx=10,winIndex=0;
	int Sy=(Sx*Ymax)/Xmax;
	float x,y;
	
	Point3f Vview(2,2,15),Vup(0,1,0);	// view vector
	Vview.Normalize();
	Point3f n2=Vview,n0=n2^Vup;
	n0.Normalize();
	Point3f n1=n0^n2;

	int d=3;
	Point3f npe=Pe,Pcenter=Pe+d*n2;
	npe.Normalize();//print(npe);
	Point3f P00=Pcenter-(Sx/2)*n0-(Sy/2)*n1,Pp;
	//vector<Sphere> allSpheres;
	//allSpheres.pb(sphere1);
	
	Sphere sphere1(Point3f(2,2,25),15, Color(1,0,0),0);
	Sphere sphere2(Point3f(2,2,5),2, Color(1,0.2,0.7),1);
	Sphere sphere3(Point3f(-1,5,8),3, Color(0.1,0.5,1),2);
	//Sphere sphere4(Point3f(1,1,6),3, Color(0,0.5,1),3);
	
	Plane plane1(Point3f(0,-1,0), Point3f(0,100,0),Color(0.2,0.2,0.2));
	Plane plane2(Point3f(-1,0,0), Point3f(100,0,0),Color(0.2,0.2,0.2));
	Plane plane3(Point3f(0,0,-1), Point3f(0,0,100),Color(0.2,0.2,0.2));
	vector<Object*> allObjects;
	allObjects.push_back(dynamic_cast<Object*>(&sphere2));
	allObjects.push_back(dynamic_cast<Object*>(&sphere1));
	allObjects.push_back(dynamic_cast<Object*>(&sphere3));
	//allObjects.push_back(dynamic_cast<Object*>(&sphere4));
	//
	allObjects.push_back(dynamic_cast<Object*>(&plane1));
	allObjects.push_back(dynamic_cast<Object*>(&plane2));
	allObjects.push_back(dynamic_cast<Object*>(&plane3));
	Ray myray(Pe,npe);
	int No=allObjects.size();	
	vector<Point3f > myinter;
	float rnd;
	// CHECK IF YOU'RE INSIDE ANY OF THE SPHERES ///////////////////////////
	for (int i = 0; i < allObjects.size(); i++)
	{
		if(!allObjects[i]->isEyeOutside(Pe))
		{	cout<<"Eye is inside the object \n";return 1;}
	}
	cout<<"Eye is outside all the object \n";
	
	// WRITE TO IMAGE //////////////////
	
	vector <Point3f> tmp;
	//tmp.resize(2);
	Color finalColor;
	Point3f interPoint(0,0,0);
	int spNum=0;
	for (int I = 0; I < Xmax; I++)
	{
		for (int J = 0; J < Ymax; J++)
		{
			index=J*Xmax + I;
			// shooting rays from center of the pizel

			x=(I+0.5)/Xmax;
			y=(J+0.5)/Ymax;

			Pp=P00+(x*Sx)*n0+(y*Sy)*n1;	//Direction to shoot the ray
			myray=Ray(Pe, Pp);	// This is the ray that we will shoot from camera to find out the color at pixel x,y
			tmp.clear();
			for (int i = 0; i < allObjects.size(); i++)
			{
				interPoint=allObjects[i]->getIntersectionPoint(myray);
				//myinter.pb(interPoint);
				//if(interPoint.x!=-1)	
				tmp.pb(interPoint);	//get intersection points
			}
			//cout<<"tmpsize: "<<tmp.size();cout<<endl;
			//if(interPoint.x == -1)

			//find which point wins
			pixels[index].r=0;
			pixels[index].g=0;
			pixels[index].b=0;

			winIndex= findWinningPointIndex(tmp , Pe);
			if(winIndex==-1)
				continue;
			//assign that index object's color to the image
			//shader equatoins:
			

			if(allObjects[winIndex]->objectName=="Sphere")	
			{
				//*
				spNum = ((Sphere*)allObjects[winIndex])->id;
				if(spNum==0)	finalColor = ((Sphere*)allObjects[winIndex])->goochShader(myray);
				else if(spNum==1)	finalColor = ((Sphere*)allObjects[winIndex])->phongShader(myray);
				else if(spNum==2)	finalColor = ((Sphere*)allObjects[winIndex])->lambertShader(myray);
				else 	finalColor = ((Sphere*)allObjects[winIndex])->lambertShader(myray);
				// */
				//finalColor = ((Sphere*)allObjects[winIndex])->phongShader(myray);

				spNum = (spNum + 1)%3;
				//cout<<spNum<<" ";
			}
			if(allObjects[winIndex]->objectName=="Plane")	
			{
				finalColor = ((Plane*)allObjects[winIndex])->lambertShader(myray);
			}
			pixels[index].r=finalColor.red;
			pixels[index].g=finalColor.green;
			pixels[index].b=finalColor.blue;
		}
	}
	time_t newTime;
	time(&newTime);
	std::string number;std::stringstream strstream;strstream << newTime;strstream >> number;
	string fileName = "scene_" + number + ".bmp";cout<<fileName;
	savebmp(fileName.c_str(),Xmax,Ymax,dpi,pixels);
	////////////////////////////////////////////////////////////////////////
	
	return 0;
	
}

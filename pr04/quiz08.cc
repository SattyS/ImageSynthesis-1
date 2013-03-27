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
/*************************************** OBJECT CLASS **********************************************************************/
class Object 
{
	public:
		Object(){}
		virtual Point3f getIntersectionPoints(Ray ray){return Point3f(0,0,0);}
		virtual bool isEyeOutside(Point3f Pe){return false;}
		virtual Point3f getColor(){return Point3f(0,0,0);}
};
/************************************** SPHERE CLASS ***********************************************************************/
class Sphere : public Object
{
	public:
	double radius;
	Point3f center,color;
	Sphere()
	{
		radius= 0.0;
		center=Point3f(0,0,0);
		color=Point3f(0,0,0);
	}
	Sphere(Point3f c, double rad,Point3f col)
	{
		radius= rad;
		center=c;
		color=col;
	}
	Sphere(Point3f c, double rad)
	{
		radius= rad;
		center=c;
		color=Point3f(1,1,1);
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
	virtual Point3f getColor(){	return color;}
	virtual Point3f getIntersectionPoints(Ray ray)
	{
		Point3f ret;
		//cout<<"assuming that eye is outside the sphere \n";
		
		double c=(center-ray.origin)%(center-ray.origin)-(radius*radius),b=(center-ray.origin)%ray.direction;
		double delta=b*b-c,t1,t2,tmin;
		if(delta<0)	{
			cout<<"delta<0\n";return Point3f(-1,-1,-1);
		}
		/*
		if(b<0)	{
			cout<<"b<0\n";return Point3f(-1,-1,-1);
		}
		// */
		t1=b+sqrt(delta);t2=b-sqrt(delta);
		cout<<t1<<" , "<<t2<<endl;
		tmin=min(t1,t2);
		cout<<"min t is : "<<tmin<<endl;
		ret=Point3f((ray.origin+tmin*ray.direction));
		//print(ret);
		//return ret;
		/*
		//tmin=min(t1,t2);
		if(b==0 || t1==t2)	//if ray is a tangent to this sphere
		{
			ret.push_back(Point3f((ray.origin+t1*ray.direction)));return ret;
		}
		//(ray.direction).Normalize();
		Point3f pt1(ray.origin + t1*ray.direction),pt2(ray.origin + t2*ray.direction);
		//cout<<2*sqrt(delta)<<" intersection points are : ";print(pt1);print(pt2);
		ret.push_back(pt1);ret.push_back(pt2);
		
		cout<<(pt1-pt2).Length();cout<<endl;
		*/
		return ret;
	}
};
int findWinningPointIndex(vector<Point3f > myinter ,Point3f Pe)
{
	float len=0,minLen=INT_MAX,index=-1;
	for (int i = 0; i < myinter.size(); i++)
	{
		if(myinter[i].x == -1)	continue;
		len=(myinter[i]-Pe).Length();
		if(len < minLen)
		{
			minLen=len;
			index=i;
		}
	}
	
	return index;
}

Point3f findArea(Point3f p0, Point3f p1, Point3f p2)
{
	Point3f v1 = p1-p0,v2=p2-p0, Area;
	Area = v1^v2;
	cout<<"main area: "<<Area.Length()<<endl; print(Area);
	return Area;

}


/********************************************* MAIN ****************************************************************/
int main (int argc, char const* argv[])
{
	
	Point3f p0(41,70,-91),p1(45,80,-88),p2(48,90,-82),ph(51.02,97.4,-79.96), v1 = p1-p0,v2=p2-p0, Area;
	
	Area= findArea(p0,p1,ph);
	Area= findArea(p0,ph,p2);
	Area= findArea(p2,ph,p1);


	//Area = v1^v2;

	//cout<<"main area: "<<Area.Length()<<endl; print(Area);
	
	return 0;
	/*
	Point3f Pe(1,1,1);	//camera or eye position
	int antiAliasing=0;
	int Xmax = 500;
	int Ymax = 500;
	int n=Ymax*Xmax,dpi=72;
	RGBType *pixels= new RGBType[n];
	int index=0,M=2,N=2; 
	int Sx=10,winIndex=0;
	int Sy=(Sx*Ymax)/Xmax;
	float x,y;
	Point3f Vview(3,3,3),Vup(0,1,0);	// view vector
	Vview.Normalize();
	Point3f n2=Vview,n0=n2^Vup;
	n0.Normalize();
	
	Point3f n1=n0^n2;
	int d=1;
	Point3f npe=Pe,Pcenter=Pe+d*n2;
	npe.Normalize();//print(npe);
	Point3f P00=Pcenter-(Sx/2)*n0-(Sy/2)*n1,Pp;
	
	//cout<<"Pcenter";print(Pcenter);
	
	//vector<Sphere> allSpheres;
	//allSpheres.pb(sphere1);
	
	Sphere sphere1(Point3f(15,15,15),15, Point3f(0.5,1,0.2));
	//Sphere sphere2(Point3f(5,5,5),2.0, Point3f(1,1,.5));
	vector<Object*> allObjects;
	allObjects.push_back(dynamic_cast<Object*>(&sphere1));
	//allObjects.push_back(dynamic_cast<Object*>(&sphere2));
	
	
	Ray myray(Pe,npe);
	
	vector<Point3f > myinter;
	float rnd;
	// CHECK IF YOU'RE INSIDE ANY OF THE SPHERES ///////////////////////////
	for (int i = 0; i < allObjects.size(); i++)
	{
		if(!allObjects[i]->isEyeOutside(Pe))
		{	cout<<"Eye is inside the object \n";return 1;}
	}
	cout<<"Eye is outside all the object \n";
	
	//allObjects.pb(sphere1);
	
	// WRITE TO IMAGE //////////////////
	
	vector <Point3f> tmp;
	//tmp.resize(2);
	
	Point3f interPoint(0,0,0);
	for (int I = 0; I < Xmax; I++)
	{
		for (int J = 0; J < Ymax; J++)
		{
			
			//for (int p = 0; p < M; p++)
			//{
			//	for (int q = 0; q < N; q++)
			//	{

			index=J*Xmax + I;

			// shooting rays from center of the pizel

			if(antiAliasing!=1)
			{
				x=(I+0.5)/Xmax;
				y=(J+0.5)/Ymax;
			}
			//else
			//{
			//	rnd=(float)((float)(rand()%100)/(100.0));
			//	x=(I+(p/M)+((rnd)/M))/Xmax;
			//	y=(J+(q/N)+((rnd)/N))/Ymax;
			//}
			//cout<<;cout<<endl;

			Pp=P00+(x*Sx)*n0+(y*Sy)*n1;	//Direction to shoot the ray
			myray=Ray(Pe, Pp);	// This is the ray that we will shoot from camera to find out the color at pixel x,y
			tmp.clear();
			for (int i = 0; i < allObjects.size(); i++)
			{
				interPoint=allObjects[i]->getIntersectionPoints(myray);
				//myinter.pb(interPoint);
				//if(interPoint.x!=-1)	
				tmp.pb(interPoint);	//get intersection points
			}

			//cout<<"tmpsize: "<<tmp.size();cout<<endl;

			//if(interPoint.x == -1)

			//find which point wins
			//if(tmp.size()==1)	winIndex=0;
			winIndex= findWinningPointIndex(tmp , Pe);
			if(winIndex==-1)
			{
				if(antiAliasing!=1)
				{
					pixels[index].r=0;
					pixels[index].g=0;
					pixels[index].b=0;
					continue;
				}
				else	continue;

			}
			//assign that index object's color to the image
			if(antiAliasing==1)
			{
				float cosTheta=((Pp-((Sphere*)allObjects[winIndex])->center)%(Pp-Pe))/((Pp-Pe).Length()*(Pp-((Sphere*)allObjects[winIndex])->center).Length());
				Point3f c = allObjects[winIndex]->getColor(), background(0,0,0),final;
				final = (cosTheta)*background + (1-cosTheta)*c;

				pixels[index].r=final.x;
				pixels[index].g=final.y;
				pixels[index].b=final.z;

			}
			else
			{
				pixels[index].r+=(allObjects[winIndex]->getColor().x)/(M*N);
				pixels[index].g+=(allObjects[winIndex]->getColor().y)/(M*N);
				pixels[index].b+=(allObjects[winIndex]->getColor().z)/(M*N);
			}
		//}
		//}
		}
	}
	savebmp("scene.bmp",Xmax,Ymax,dpi,pixels);
	////////////////////////////////////////////////////////////////////////
	
	return 0;
*/	
}

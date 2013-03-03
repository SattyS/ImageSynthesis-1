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
#include "ObjectsDefine.h"

using namespace cy;
using namespace std;
#define print(p) printf("(%f,%f,%f) \n",p.x,p.y,p.z);
#define pb(a) push_back(a)
Point3f Pe(0,0,0);      //camera or eye position
Point3f PL(0,20,10);
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
	
	Point3f Vview(0,0,1),Vup(0,1,0);	// point the view vector to focus on a particular point from Pe
	//Vview = Vview - Pe;
	Vview.Normalize();
	//Vup.Normalize();
	Point3f n2=Vview,n0=n2^Vup;
	n0.Normalize();
	Point3f n1=n0^n2;
	n1.Normalize();
	int d=3;
	Point3f npe=Pe,Pcenter=Pe+d*n2;
	npe.Normalize();//print(npe);
	Point3f P00=Pcenter-(Sx/2)*n0-(Sy/2)*n1,Pp;
	//vector<Sphere> allSpheres;
	//allSpheres.pb(sphere1);
	
	Sphere sphere1(Point3f(0,-5,20),10, Color(1,0.8,0),1);
	//Sphere sphere2(Point3f(2,2,5),2, Color(1,0.2,0.7),1);
	//Sphere sphere3(Point3f(-1,5,8),3, Color(0.1,0.5,1),2);
	//Sphere sphere4(Point3f(1,1,6),3, Color(0,0.5,1),3);
	
	Plane plane1(Point3f(0,-1,0), Point3f(0,100,0), Color(1,0,1));
	Plane plane2(Point3f(-1,0,0), Point3f(100,0,0), Color(0.2,1,0.2));
	Plane plane3(Point3f(0,0,-1), Point3f(0,0,100), Color(1,0.2,0.2));
	vector<Object*> allObjects;
	allObjects.push_back(dynamic_cast<Object*>(&sphere1));
	//allObjects.push_back(dynamic_cast<Object*>(&sphere2));
	//allObjects.push_back(dynamic_cast<Object*>(&sphere3));
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
			Point3f shadowInter,rayStart = tmp[winIndex];
			int flag=0;
			Ray jujuRay(rayStart , PL-rayStart);
			flag=0;
			tmp.clear();

			for(int j=0;j<allObjects.size();j++)
			{
				if(j==winIndex)	continue;
				shadowInter = allObjects[j]->getIntersectionPoint(jujuRay);
				tmp.pb(shadowInter);

			}
			int mywinIndex= findWinningPointIndex(tmp , PL);
			if(mywinIndex !=-1 && tmp[mywinIndex]==rayStart)
				continue;
			else if((mywinIndex != -1 && (rayStart-tmp[mywinIndex]).Length() > (rayStart - PL).Length()) || (rayStart-tmp[mywinIndex])%(rayStart - PL)<0 )
				continue;
			
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
				//finalColor = ((Plane*)allObjects[winIndex])->getColor();
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

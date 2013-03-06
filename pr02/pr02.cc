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
class SpotLight
{
    public:
      Point3f source,direction;
      double angle;

      SpotLight(){  source=Point3f(0,0,0),direction=Point3f(0,0,0);angle = 30;}
      SpotLight(Point3f src,Point3f dir,double an){  source=src,direction=dir;angle = an;}
      SpotLight(Point3f src,Point3f dir,float an){  source=src,direction=dir;angle = an;}


};

Point3f Pe(0,0,0);      //camera or eye position
SpotLight spotLight(Point3f(0,14,25),Point3f(0,0,1),160.0/180.0);

//Point3f PL = spotLight.source;
Point3f DirectionLight(0,-1,0);

/********************************************* MAIN ****************************************************************/
int main (int argc, char const* argv[])
{
	int antiAliasing=0;
	
        bool spotlightEnabled=false;
	
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
	
	Sphere sphere1(Point3f(0,0,25),10, Color(1,0.8,0),1);
	Sphere sphere2(Point3f(-3,-2,10),3, Color(1,0.2,0.7),2);
	//Sphere sphere3(Point3f(-1,5,8),3, Color(0.1,0.5,1),3);
	//Sphere sphere4(Point3f(1,1,6),3, Color(0,0.5,1),4);
	
	Plane plane1(Point3f(0,-1,0), Point3f(0,30,0), Color(.3,0,.3));
	Plane plane2(Point3f(-1,0,0), Point3f(30,0,0), Color(.3,1,.3));
	Plane plane3(Point3f(0,0,-1), Point3f(0,0,60), Color(.3,1,1));
	Plane plane4(Point3f(0,0,1), Point3f(0,0,-30), Color(0,1,0));
	Plane plane5(Point3f(0,1,0), Point3f(0,-30,0), Color(1,0,0));
	Plane plane6(Point3f(1,0,0), Point3f(-30,0,0), Color(1,1,0));

	vector<Object*> allObjects;
	//allObjects.push_back(dynamic_cast<Object*>(&sphere1));
	allObjects.push_back(dynamic_cast<Object*>(&sphere2));
	//allObjects.push_back(dynamic_cast<Object*>(&sphere3));
	//allObjects.push_back(dynamic_cast<Object*>(&sphere4));
	//
	allObjects.push_back(dynamic_cast<Object*>(&plane1));
	allObjects.push_back(dynamic_cast<Object*>(&plane2));
	allObjects.push_back(dynamic_cast<Object*>(&plane3));
	allObjects.push_back(dynamic_cast<Object*>(&plane4));
	allObjects.push_back(dynamic_cast<Object*>(&plane5));
	allObjects.push_back(dynamic_cast<Object*>(&plane6));
	Ray myray(Pe,npe);
	int No=allObjects.size();	
	vector<Point3f > myinter;
	float rnd;
        Point3f PL;
	// CHECK IF YOU'RE INSIDE ANY OF THE SPHERES ///////////////////////////
	for (int i = 0; i < allObjects.size(); i++)
	{	
          if(!allObjects[i]->isEyeOutside(Pe) )
		{	cout<<"Eye is inside the object \n";return 1;}
	}
	cout<<"EEye and light are outside all the object \n";
	
        vector<Point3f> lights;
        lights.pb(Point3f(10,24,30));
        lights.pb(Point3f(-10,24,30));
        lights.pb(Point3f(0,24,0));

        double alpha0= cos(spotLight.angle);
        int softShadowFlag=0;
	vector <Point3f> tmp;
	//tmp.resize(2);
	Color finalColor(0,0,0);
        Color TotalfinalColor(0,0,0);
        double intensity=1.0/(lights.size());
	Point3f interPoint(0,0,0);
	int spNum=0;
        double alpha=0;
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
				//if((allObjects[winIndex]->objectName=="Plane"))
				//	cout<<"intersection point with plane : "; print(interPoint);
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
                        Point3f NH ;
                        finalColor = Color(0,0,0);



                        

                       TotalfinalColor = Color(0,0,0);


                        int numLights = lights.size();
                        //cout<<numLights<<endl;
                        for(int lightNum = 0; lightNum < numLights; lightNum++)
                        {
		                  PL = lights[lightNum];

		                  //print(PL);
		                  finalColor = Color(0,0,0);
		                  Ray jujuRay(rayStart , PL-rayStart);

		                  flag=0;
		                  tmp.clear();
		                  double shadowDist=0.0,distToPL=(PL-rayStart).Length(),ratio=1.0,o=1,maxDark=1;
		                  int power=1;
		                  double cDirect;
		                  //*/ FOR SHARP SHADOWS 
		                  // * get all other intersection points with jujuRay except for the object itself (winIndex)
		                  for(int j=0;j<allObjects.size();j++)
		                  {
		                    if(j==winIndex)	continue;
		                    // for soft shadows, find the distance between the intersection points also
		                    if(allObjects[j]->objectName=="Sphere" )	
		                    {
		                      shadowDist += (((Sphere*)allObjects[j])->getTwoDelta(jujuRay))/(2*( ((Sphere*)allObjects[j])->radius)) ;
		                      shadowDist = min(maxDark, shadowDist);
		                    }
		                    shadowInter = allObjects[j]->getIntersectionPoint(jujuRay);
		                    tmp.pb(shadowInter);

		                  }
		                  int mywinIndex= findWinningPointIndex(tmp , rayStart);
		                  if(mywinIndex !=-1 && allObjects[mywinIndex]->objectName=="Sphere")	
		                  {
		                    NH  = (tmp[mywinIndex] - ((Sphere*)(allObjects[mywinIndex]))->center   );
		                  }
		                  //*
		                  else if(mywinIndex !=-1 && allObjects[mywinIndex]->objectName=="Plane")	
		                  {
		                    NH = (-1)*(((Plane*)(allObjects[mywinIndex]))->normalVector);
		                  }


		                  NH.Normalize();
		                  DirectionLight.Normalize();
		                  cDirect = DirectionLight%NH;                         // for directional light control value
		                  //cout<<cDirect<<endl;
		                  if(cDirect < 0.0 )
		                  {
		                    //cout<<cDirect<<endl;
		                    cDirect = 1;
		                  }
		                  else {
		                    //cout<<cDirect<<endl;
		                    cDirect = 0;
		                  }

		                  // */
		                  // uncomment the following line for normal point light
		                  cDirect = 0;
		                  if(mywinIndex !=-1)// && tmp[mywinIndex]!=rayStart)
		                  {
		                    // if distance between startRay and the point tmp[mywinIndex] < dist b/n startRay and PL then black is the color -> continue
		                    double pointToInter = (rayStart-tmp[mywinIndex]).Length(), pointToLight = ( rayStart-PL).Length();

		                    if( pointToInter < pointToLight)
		                    {
		                      //ratio =  1 - (shadowDist/distToPL)  ;

		                      ratio =  1 - pow(shadowDist,power)  ;
		                      // comment the following line to have sharp shadows
		                      softShadowFlag = 1;


		                    }
		                  }
		                  // FOR Soft shadows
		                  Color shadowColor(ratio,ratio,ratio);
		                  //shadowColor=shadowColor*ratio;


		                  alpha = ((rayStart - spotLight.source)%spotLight.direction)/((rayStart - spotLight.source).Length()*(spotLight.direction).Length());
		                  if(alpha >alpha0) alpha=1;
		                  else alpha = 0;

		                  if(allObjects[winIndex]->objectName=="Sphere")	
		                  {
		                    if(spotlightEnabled)
		                      finalColor = finalColor+ (((Sphere*)allObjects[winIndex])->phongShader(myray,PL))*alpha;
		                    else
		                      finalColor = finalColor+ ((Sphere*)allObjects[winIndex])->phongShader(myray,PL);

		                    if(mywinIndex!=-1 && softShadowFlag==0)
		                      finalColor = finalColor+  Color(0,0,0);

		                    if(cDirect==0)
		                      finalColor = finalColor*shadowColor;

		                                             }
		                  if(allObjects[winIndex]->objectName=="Plane")	
		                  {
		                    if(spotlightEnabled)
		                      finalColor = finalColor+ (((Plane*)allObjects[winIndex])->phongShader(myray,PL))*alpha;
		                    else
		                    {
		                      //finalColor = finalColor+ ((Plane*)allObjects[winIndex])->lambertShader(myray,PL);

		                      finalColor = finalColor+  ((Plane*)allObjects[winIndex])->getColor();
		                    }


		                    //finalColor = finalColor+  ((Plane*)allObjects[winIndex])->getColor();
		                    if(mywinIndex!=-1 && softShadowFlag==0)
		                      finalColor = finalColor+  Color(0,0,0);
		                    if(cDirect==0)
		                      finalColor =  finalColor*shadowColor;
		                    
		                    //printf("Final color in iteration %i : \n", lightNum);
		                    //finalColor.printColor();


		                  }

		                  TotalfinalColor = TotalfinalColor + finalColor*intensity;
                                  //cout<<"total color: "; TotalfinalColor.printColor();cout<<endl;
                          

                        }
			finalColor = TotalfinalColor;
                        /*
                        if(finalColor.red > 1.0)  finalColor.red = 1.0;
                        if(finalColor.green > 1.0)  finalColor.green = 1.0;
                        if(finalColor.blue > 1.0)  finalColor.blue = 1.0;

                        if(finalColor.red < 0)  finalColor.red = 0;
                        if(finalColor.green < 0)  finalColor.green = 0;
                        if(finalColor.blue < 0)  finalColor.blue = 0;
                        // */
			

			
                        

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

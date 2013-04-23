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

#include <string>
#include <vector>
#include <sstream>
#include <fstream>

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
#include<limits.h>
#include "ObjectsDefine.h"
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
SpotLight spotLight(Point3f(0,10,5),Point3f(0,0,1),60.0/180.0);
//Enable each variable to enable textures on them
bool sphereTextureEnabled=false, genericTextureEnabled=false,planeTextureEnabled=true, textureRefractionMapEnabled=false;
bool glossyEnabled=0,translucencyEnabled=1;
//Point3f PL = spotLight.source;
Point3f DirectionLight(0,-1,0);

/*** READ PPM ******/
/********************************************* MAIN ****************************************************************/
/*
void printVector(obj_vector *v)
{
	printf("%.2f,", v->e[0] );
	printf("%.2f,", v->e[1] );
	printf("%.2f ", v->e[2] );
}

*/
int numRecursion = 4;
Color rayTracer(Ray myray, Point3f PL, vector<Object*> allObjects , int depth, bool spotlightEnabled ,int softShadowFlag);
int main (int argc, char const* argv[])
{
	int antiAliasing=1;
    bool spotlightEnabled=false;
	
	int Xmax = 500;
	int Ymax = 500;
	int n=Ymax*Xmax,dpi=72;

	readPPM();
	readPPM1();
	readPPM2();
	readPPMBumpMap();
	
	RGBType *pixels= new RGBType[n];
	int index=0,M=2,N=M; 
	int Sx=10,winIndex=0;
	int Sy=(Sx*Ymax)/Xmax;
	float x,y;
	
	Point3f Vview(0,0,65),Vup(0,1,0);	// point the view vector to focus on a particular point from Pe
	//Vview = Point3f(0,0,1);
	//Vview = Vview - Pe;
	Vview.Normalize();
	//Vup.Normalize();
	Point3f n2=Vview,n0=n2^Vup;
	n0.Normalize();
	Point3f n1=n0^n2;
	n1.Normalize();
	int d=10;
	Point3f npe=Pe,Pcenter=Pe+d*n2;
	npe.Normalize();//print(npe);
	Point3f P00=Pcenter-(Sx/2)*n0-(Sy/2)*n1,Pp;

	//vector<Sphere> allSpheres;
	//allSpheres.pb(sphere1);
	//string genericObjFileName="cube1.obj";
	//printf("%s", genericObjFileName.c_str());
	
	string genericObjFileName="cube_00.obj";
	genericObjFileName = "cube_oriented.obj";
        //char *objfilename = "tetrahedron.obj";
        Sphere sphere1(Point3f(-15,0,65),15, Color(0,1,1),1,0 , 2);
	Sphere sphere2(Point3f(15,0,110),17, Color(1,0.2,0.3),2,0, 1);
	Sphere sphere3(Point3f(0,0,85),16, Color(0.1,0.5,1),3,0, 1);
	//Sphere sphere4(Point3f(1,1,6),3, Color(0,0.5,1),4);
	
	Plane plane1(Point3f(0,-1,0), Point3f(0,40,0), Color(1,1,1), "roof", 0,0);
	Plane plane2(Point3f(-1,0,0), Point3f(60,0,0), Color(1,0,0), "left", 0,0);
	Plane plane3(Point3f(0,0,-1), Point3f(0,0,60), Color(1,1,1), "front",0,0);
	Plane plane4(Point3f(0,0,1), Point3f(0,0,-80), Color(0,0,0), "back", 0,0);
	Plane plane5(Point3f(0,1,0), Point3f(0,-50,0), Color(1,1,1), "floor",0,0);
	Plane plane6(Point3f(1,0,0), Point3f(-60,0,0), Color(0,1,0), "right",0,0);

	//cout<<"debug/////";
	GenericObject cube(genericObjFileName,0.2,1.33);

	//ObjMesh triMesh= LoadObjMesh(genericObjFileName);

	//GenericObject cube1(triMesh);
	//GenericObject cube1(genericObjFileName);
	//cout<<"debug/////";

	vector<Object*> allObjects;
	//allObjects.push_back(dynamic_cast<Object*>(&cube));
	
	allObjects.push_back(dynamic_cast<Object*>(&sphere1));
	allObjects.push_back(dynamic_cast<Object*>(&sphere2));
	allObjects.push_back(dynamic_cast<Object*>(&sphere3));
	//allObjects.push_back(dynamic_cast<Object*>(&sphere4));
	
	//allObjects.push_back(dynamic_cast<Object*>(&plane1));
	//allObjects.push_back(dynamic_cast<Object*>(&plane2));
	//allObjects.push_back(dynamic_cast<Object*>(&plane3));
	//allObjects.push_back(dynamic_cast<Object*>(&plane4));
	//allObjects.push_back(dynamic_cast<Object*>(&plane5));
	//allObjects.push_back(dynamic_cast<Object*>(&plane6));

	//*/
	Ray myray(Pe,npe);
	int No=allObjects.size();	
	vector<Point3f > myinter;
	float rnd;

			//printf("psi: %f 
    Point3f PL;
	PL = spotLight.source;

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
	bool isOnePicture =true;
        for (float ior = 0; ior<=0;)
        {
          //allObjects[0]->eta = ior;
            for (int I = 0; I < Xmax; I++)
	{
		for (int J = 0; J < Ymax; J++)
		{
			// shooting rays from center of the pizel
			for (int p = 0; p < M; p++)
			{
				for (int q = 0; q < N; q++)
				{
			x=(I+0.5)/Xmax;
			y=(J+0.5)/Ymax;
			if(antiAliasing!=1)
			{
				x=(I+0.5)/Xmax;
				y=(J+0.5)/Ymax;
			}
			else
			{
				rnd=(float)((float)(rand()%100)/(100.0));
				x=(I+(p/M)+((rnd)/M))/Xmax;
				y=(J+(q/N)+((rnd)/N))/Ymax;
			}
			index=J*Xmax + I;
			Pp=P00+(x*Sx)*n0+(y*Sy)*n1;	//Direction to shoot the ray
			myray=Ray(Pe, Pp);	// This is the ray that we will shoot from camera to find out the color at pixel x,y
			
			finalColor = rayTracer(myray, PL,  allObjects , 0,spotlightEnabled ,softShadowFlag);
                        
                        Point3f refLectedRayDirection =myray.direction ;refLectedRayDirection.Normalize();
                        Color colorFromReflectedObject = finalColor;
			// comment the following line for ENVIRONMENT MAP:

                        //*                         
                        if(colorFromReflectedObject.red ==0 && colorFromReflectedObject.green ==0 && colorFromReflectedObject.blue==0)
                        {
                          double s0=1;
                          double X = refLectedRayDirection.x/s0 , Y = refLectedRayDirection.y/s0, Z = refLectedRayDirection.z/s0;

                          double psi = acos(Z);
                          //double theta = acos( Y/(float)(sqrt((1-(Z*Z))) )  );
			  double theta;// = acos((double)( Y/(double)(sqrt((1.0-(Z*Z))) )  ) );

			  if(( Y/( sqrt((1-(Z*Z)) ))>1.0 ) ||  (Y/( sqrt((1-(Z*Z)) ))<-1.0 )) 
			  {
				  //cout<<"something is happenning============================================================================ \n";
				  if(( Y/( sqrt((1-(Z*Z)) ))>1.0 )){
					  theta = acos(1.0);
				  }
				  else if(( Y/( sqrt((1-(Z*Z)) ))<-1.0 )){
					  theta = acos(-1.0);
				  }
				  //theta = asin(X/( sqrt((1-(Z*Z))) )  );
			  }
			  else
				  theta = acos((double)( Y/(double)(sqrt((1.0-(Z*Z))) )  ) );

			  double PI = 3.14;
			  double v = psi/PI, u = theta/(2*PI);

                          if(X<0)	{u = 1-u;//v=1-v;//cout<<"adfioubnwirgnw";
                          }
                          //if(v<0)	v = v+1;

                          //if(( (X>0 && X<1) && (Y>0 && Y<1) ) )   {                                  u = X; v = Y;
                          u = u*projectionImageWidth,v=v*projectionImageHeight;
                          int pixmapIndex = abs((int)v * projectionImageWidth + (int)u) * 3;

                          // printf("psi: %f , theta: %f , u: %f , v: %f \n",psi, theta, u , v);
                          //cout<<(int)( (Y * projectionImageWidth + X) * 3 )<<endl; 
                          //cout<< (float)pixmap[pixmapIndex]<<endl;
                          colorFromReflectedObject.red = (float)(pixmap[pixmapIndex])/maxcolor;
                          colorFromReflectedObject.green =(float)(pixmap[pixmapIndex+1])/maxcolor;
                          colorFromReflectedObject.blue = (float)(pixmap[pixmapIndex + 2])/maxcolor;

                          finalColor = colorFromReflectedObject;


                        }	
                        // */
			if(antiAliasing==1)
			{
				pixels[index].r=finalColor.red;  
				pixels[index].g=finalColor.green;
				pixels[index].b=finalColor.blue; 

			}
			else
			{
				pixels[index].r+=finalColor.red/(M*N);
				pixels[index].g+=finalColor.green/(M*N);
				pixels[index].b+=finalColor.blue/(M*N);
			}
			}
			}
			//pixels[index].r=finalColor.red;
			//pixels[index].g=finalColor.green;
			//pixels[index].b=finalColor.blue;
		}
	}
        
	time_t newTime;
	time(&newTime);

	std::string number;std::stringstream strstream;strstream <<newTime;strstream >> number;
	string fileName = "scene_" + number + ".bmp";cout<<fileName<<endl;
	savebmp(fileName.c_str(),Xmax,Ymax,dpi,pixels);
	////////////////////////////////////////////////////////////////////////
        ior+=0.05;
        }
	return 0;
	
}

Color rayTracer(Ray myray, Point3f PL, vector<Object*> allObjects , int depth,  bool spotlightEnabled=false ,int softShadowFlag=0)
{
//	cout<<"depth: "<<depth<<endl;
    Point3f interPoint(0,0,0);
    vector <Point3f> tmp;
	Color finalColor(0,0,0);
	Color TotalfinalColor(0,0,0);
    //double intensity=1.0/(lights.size());
	int spNum=0;
    double alpha=0;
    tmp.clear();
    double alpha0= cos(spotLight.angle);
    int winIndex=0;
    for (int i = 0; i < allObjects.size(); i++)
    {
	    interPoint=allObjects[i]->getIntersectionPoint(myray);
	    //if((allObjects[winIndex]->objectName=="genericObject"))
	    //	cout<<"intersection point with plane : "; print(interPoint);
	    //myinter.pb(interPoint);
	    //if(interPoint.x!=-1)	
	    tmp.pb(interPoint);	//get intersection points
    }
    //cout<<"tmpsize: "<<tmp.size();cout<<endl;
    //if(interPoint.x == -1)

    //find which point wins

    winIndex= findWinningPointIndex(tmp , Pe);
    if(winIndex==-1)
	    return Color(0,0,0);
    //assign that index object's color to the image
    //shader equatoins:
    Point3f shadowInter,rayStart = tmp[winIndex], interSectionPoint = tmp[winIndex];
    int flag=0;
    Point3f NH ;
    finalColor = Color(0,0,0);

    Color colorFromReflectedObject(0,0,0);
    Color colorFromRefractedObject(0,0,0);
    double ks=allObjects[winIndex]->KS;
    float rnd1 =-1 + 2*(rand()/float(RAND_MAX)),rnd2 = -1 + 2*(rand()/float(RAND_MAX)),rnd3 = -1 + 2*(rand()/float(RAND_MAX)) ;
    //rnd1 = pow(rnd1,4);
    //rnd2 = pow(rnd2,4);
    //rnd3 = pow(rnd3,4);
    //printf(" rnd1,2,3: %f, %f, %f\n ",rnd1,rnd2,rnd3);

    
    if(depth<numRecursion && ks!=0)
    {
	    Point3f n=allObjects[winIndex]->getNormal(interSectionPoint);n.Normalize();
	    Point3f v= interSectionPoint-Pe;v.Normalize();
	    v = myray.direction;v.Normalize();v=-1*v;
            // for reflection
	    // for glossy 
	    float glossyS = 0;
	    if(allObjects[winIndex]->objectName=="Sphere" )	
	    {
		    if(((Sphere*)allObjects[winIndex])->id==2)
			    glossyS=0.5;
		    if(((Sphere*)allObjects[winIndex])->id==3)
			    glossyS=0.02;
	    }
	    Point3f Vrand(rnd1,0,0);
            //if(rnd1 >1.0 || rnd1 <-1.0) cout<<"wrond rnd1!!!!!!! "<<rnd1<<endl;
	    Vrand = glossyS*Vrand;
	    Point3f refLectedRayDirection =-1*v + 2*(n%v)*n;
	    refLectedRayDirection.Normalize();
	    Ray refLectedRay(interSectionPoint, refLectedRayDirection);
            //colorFromReflectedObject = rayTracer(refLectedRay, PL,  allObjects , depth+1, spotlightEnabled ,softShadowFlag); 
	    if(glossyEnabled)
	    {
              float NN=2;
              for(int i=1;i<NN;i++)
              {
                rnd1 =-1 + 2*(rand()/float(RAND_MAX));
                Vrand = Point3f(rnd1,0,0);

                Vrand = glossyS*Vrand;
		    refLectedRayDirection = refLectedRayDirection + Vrand;
		    refLectedRayDirection.Normalize();
		    refLectedRay = Ray(interSectionPoint, refLectedRayDirection);
		    colorFromReflectedObject = (colorFromReflectedObject + rayTracer(refLectedRay, PL,  allObjects , depth+1, spotlightEnabled ,softShadowFlag))/NN; 
              }
	    }
	    else
		    colorFromReflectedObject = (colorFromReflectedObject + rayTracer(refLectedRay, PL,  allObjects , depth+1, spotlightEnabled ,softShadowFlag))/1.0; 

	    // comment the following line for ENVIRONMENT MAP:
	    //*
              if(colorFromReflectedObject.red ==0 && colorFromReflectedObject.green ==0 && colorFromReflectedObject.blue==0)
	    {
		    double s0=1;
                    double X = refLectedRayDirection.x/s0 , Y = refLectedRayDirection.y/s0, Z = refLectedRayDirection.z/s0;

                    double psi = acos(Z);
		    double theta;// = acos((double)( Y/(double)(sqrt((1.0-(Z*Z))) )  ) );

		    if(( Y/( sqrt((1-(Z*Z)) ))>1.0 ) ||  (Y/( sqrt((1-(Z*Z)) ))<-1.0 )) 
		    {
			    //cout<<"something is happenning============================================================================ \n";
			    if(( Y/( sqrt((1-(Z*Z)) ))>1.0 )){
				    theta = acos(1.0);
			    }
			    else if(( Y/( sqrt((1-(Z*Z)) ))<-1.0 )){
				    theta = acos(-1.0);
			    } 			    //theta = asin(X/( sqrt((1-(Z*Z))) )  );
		    }
		    else
			    theta = acos((double)( Y/(double)(sqrt((1.0-(Z*Z))) )  ) );
                    double PI = 3.14;
                    double v = psi/PI, u = theta/(2*PI);

                    if(X<0)	{u = 1-u;//v=1-v;//cout<<"adfioubnwirgnw";
                    }
		    //if(v<0)	v = v+1;

		    //if(( (X>0 && X<1) && (Y>0 && Y<1) ) )   {                                  u = X; v = Y;
		    u = u*projectionImageWidth,v=v*projectionImageHeight;
		    int pixmapIndex = abs((int)v * projectionImageWidth + (int)u) * 3;

		    //printf("psi: %f , theta: %f , u: %f , v: %f \n",psi, theta, u , v);
		    //printf("X: %f , Y: %f , Z: %f, theta: %e \n",X, Y, Z,( Y/(sqrt((1-(Z*Z))) )  ));
		    //cout<<(int)( (Y * projectionImageWidth + X) * 3 )<<endl; 
		    //cout<< (float)pixmap[pixmapIndex]<<endl;
		    colorFromReflectedObject.red = (float)(pixmap[pixmapIndex])/maxcolor;
		    colorFromReflectedObject.green =(float)(pixmap[pixmapIndex + 1])/maxcolor;
		    colorFromReflectedObject.blue = (float)(pixmap[pixmapIndex + 2])/maxcolor;
            }
            // */


    }
    double eta;
    if(textureRefractionMapEnabled)
      eta = allObjects[winIndex]->findEta(interSectionPoint);
    else
      eta = allObjects[winIndex]->eta;

    if(depth<numRecursion && eta>0.05)
    {
	    Point3f n=allObjects[winIndex]->getNormal(interSectionPoint);n.Normalize();
	    Point3f v= interSectionPoint-Pe;v.Normalize();
	    v = myray.direction;v.Normalize();v=-1*v;
            // for refraction
            float c = n%v;
            float a = -1.0/eta;
            float b = (c-sqrt(c*c - 1 + eta*eta))/eta;

            //printf("eta: %f, a: %f , b: %f, c: %f\n" ,eta,a,b,c );

            Point3f refractedRayDirection =-1*a*v -(1-a)*n;refractedRayDirection.Normalize();
            Ray refractedRay(interSectionPoint, refractedRayDirection);

            //cout<<"interSectionPoint: "; print(interSectionPoint);
	    Point3f nextIntersection = allObjects[winIndex]->getOtherIntersectionPoint(refractedRay);
            Point3f nnew=allObjects[winIndex]->getNormal(nextIntersection);nnew.Normalize();
            Point3f vnew= interSectionPoint-Pe;vnew.Normalize();
	    vnew = refractedRay.direction;vnew.Normalize();vnew=-1*vnew;
            // for refraction
            float etaout=1.0/eta;
            float cnew = nnew%vnew;
            float anew = -1.0/etaout;
            float bnew = (cnew-sqrt(cnew*cnew - 1 + etaout*etaout))/etaout;

            //printf("new: eta: %f, a: %f , b: %f, c: %f\n" ,eta,a,b,c );
            Point3f refractedRayDirectionnew =-1*(1-anew)*vnew - anew*nnew;refractedRayDirectionnew.Normalize();
            Ray refractedRaynew(nextIntersection, refractedRayDirectionnew);
            if(nextIntersection.x ==-1 && nextIntersection.y ==-1 && nextIntersection.z ==-1 )
            {
              cout<<"Were dooomed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
              cout<<"interSectionPoint: "; print(interSectionPoint);
            }
            //cout<<"nextIntersection:";print(nextIntersection);
            //cout<<"refractedRaynew : "; print(refractedRaynew.direction);
            // */
            // for translucency 
            float transS = 0;
	    if(allObjects[winIndex]->objectName=="Sphere" )	
	    {
		    if(((Sphere*)allObjects[winIndex])->id==2)
			    transS=0;
		    if(((Sphere*)allObjects[winIndex])->id==3)
			    transS=0.08;
	    }
	    Point3f Vrand(rnd1,0,0);
            //if(rnd1 >1.0 || rnd1 <-1.0) cout<<"wrond rnd1!!!!!!! "<<rnd1<<endl;
	    Vrand = transS*Vrand;
            //colorFromReflectedObject = rayTracer(refLectedRay, PL,  allObjects , depth+1, spotlightEnabled ,softShadowFlag); 
	    if(translucencyEnabled)
	    {
              float NN=2;
              for(int i=1;i<NN;i++)
              {
                rnd1 =-3 + 6*(rand()/float(RAND_MAX));
                Vrand = Point3f(rnd1,rnd2,rnd3);

                Vrand = transS*Vrand;
                refractedRayDirectionnew = refractedRayDirectionnew + Vrand;
                refractedRayDirectionnew.Normalize();
                refractedRaynew = Ray(nextIntersection, refractedRayDirectionnew);
                colorFromRefractedObject= (colorFromRefractedObject + rayTracer(refractedRaynew, PL,  allObjects , depth+1, spotlightEnabled ,softShadowFlag))/NN; 
              }
	    }
	    else
              colorFromRefractedObject = rayTracer(refractedRaynew, PL,  allObjects , depth+1, spotlightEnabled ,softShadowFlag);
	
            //cout<<"colorFromRefractedObject: ";
            //colorFromRefractedObject.printColor();
            //cout<<endl; 
	    // comment the following line for ENVIRONMENT MAP:
	    //*
            if(colorFromRefractedObject.red ==0 && colorFromRefractedObject.green ==0 && colorFromRefractedObject.blue==0)
	    {
		    double s0=1;
                    double X = refractedRayDirectionnew.x/s0 , Y = refractedRayDirectionnew.y/s0, Z = refractedRayDirectionnew.z/s0;

                    double psi = acos(Z);
                    double theta = acos((double)( Y/(double)(sqrt((1.0-(Z*Z))) )  ) );
                    if(( Y/( sqrt((1-(Z*Z)) ))>1.0 ) ||  (Y/( sqrt((1-(Z*Z)) ))<-1.0 )) 
                    {
			    //cout<<"something is happenning============================================================================ \n";
			    if(( Y/( sqrt((1-(Z*Z)) ))>1.0 )){
				    theta = acos(1.0);
			    }
			    else if(( Y/( sqrt((1-(Z*Z)) ))<-1.0 )){
				    theta = acos(-1.0);
			    } 			    //theta = asin(X/( sqrt((1-(Z*Z))) )  );
		    }
		    else
			    theta = acos((double)( Y/(double)(sqrt((1.0-(Z*Z))) )  ) );



                    double PI = 3.14;
                    double v = psi/PI, u = theta/(2*PI);

                    if(X<0)	{u = 1-u;//v=1-v;//cout<<"adfioubnwirgnw";
                    }
		    //if(v<0)	v = v+1;

		    //if(( (X>0 && X<1) && (Y>0 && Y<1) ) )   {                                  u = X; v = Y;
		    u = u*projectionImageWidth,v=v*projectionImageHeight;
		    int pixmapIndex = abs((int)v * projectionImageWidth + (int)u) * 3;

		    //printf("psi: %f , theta: %f , u: %f , v: %f \n",psi, theta, u , v);
		    //printf("X: %f , Y: %f , Z: %f, theta: %e \n",X, Y, Z,( Y/(sqrt((1-(Z*Z))) )  ));
		    //cout<<(int)( (Y * projectionImageWidth + X) * 3 )<<endl; 
		    //cout<< (float)pixmap[pixmapIndex]<<endl;
		    colorFromRefractedObject.red = (float)(pixmap[pixmapIndex])/maxcolor;
		    colorFromRefractedObject.green =(float)(pixmap[pixmapIndex + 1])/maxcolor;
		    colorFromRefractedObject.blue = (float)(pixmap[pixmapIndex + 2])/maxcolor;
            }
            // */


    }


    TotalfinalColor = Color(0,0,0);


    //int numLights = lights.size();
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
	    else if(allObjects[j]->objectName=="genericObject" )	
	    {
		    shadowDist += (((GenericObject*)allObjects[j])->twoDist(jujuRay))/5.0;
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
	    {
		    double X,Y,S0=10,S1=10,Z;
		    S0=200;
		    S1=S0;


		    if(!sphereTextureEnabled)	{
			    finalColor = finalColor + (((Sphere*)allObjects[winIndex])->phongShader(myray,PL));
			    
			    if(((Sphere*)allObjects[winIndex])->id==2)
			    {	
				    //cout<<"finalColor before: ";finalColor.printColor();cout<<endl;
			    }
			    finalColor = finalColor*shadowColor;
                            //cout<<"eta: "<<eta<<endl;
                            //cout<<"colorFromRefractedObject: "; colorFromRefractedObject.printColor();cout<<endl;
			    if(ks!=0 || eta!=0)
                            finalColor =(colorFromRefractedObject*(1-ks) + colorFromReflectedObject*ks);
                            //cout<<"finalColor: "; finalColor.printColor();cout<<endl;
			    /*
                             * float myalpha = 0.3;
                            if(depth < numRecursion && eta>0.05){
                              finalColor = finalColor*myalpha +colorFromRefractedObject*(1-myalpha);
                              finalColor = colorFromRefractedObject;
                            }
                            */


			    if(((Sphere*)allObjects[winIndex])->id==2)
			    {
				    //cout<<"finalColor: ";finalColor.printColor();cout<<endl;			   
			    }
			    return finalColor;
		    }
		    interSectionPoint =( interSectionPoint - ((Sphere*)allObjects[winIndex])->center);
		    interSectionPoint = interSectionPoint/((Sphere*)allObjects[winIndex])->radius;

		    X = interSectionPoint%Point3f(1,0,0);
		    Y = interSectionPoint%Point3f(0,1,0);
		    Z = interSectionPoint%Point3f(0,0,1);

		    /*
		       double u =X- (int)X,v = Y - (int)Y, w= Z-(int)Z;
		       if(u<0)	u = u+1;
		       if(v<0)	v = v+1;
		       if(w<0)	w = w+1;
		       */
		    double u =X- (int)X,v = Y - (int)Y;
		    if(u<0)	u = u+1;
		    if(v<0)	v = v+1;

		    //if(( (X>0 && X<1) && (Y>0 && Y<1) ) )   {                                  u = X; v = Y;
		    u = u*projectionImageWidth,v=v*projectionImageHeight;
		    int pixmapIndex = abs((int)v * projectionImageWidth + (int)u) * 3;

		    //printf("psi: %f , theta: %f , u: %f , v: %f \n",psi, theta, u , v);
		    //cout<<(int)( (Y * projectionImageWidth + X) * 3 )<<endl; 
		    //cout<< (float)pixmap[pixmapIndex]<<endl;
		    finalColor.red = (float)(pixmap[pixmapIndex])/maxcolor;
		    finalColor.green =(float)(pixmap[pixmapIndex + 1])/maxcolor;
		    finalColor.blue = (float)(pixmap[pixmapIndex + 2])/maxcolor;


		    finalColor = finalColor*(((Sphere*)allObjects[winIndex])->phongShader(myray,PL));
	    }

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
			    finalColor = finalColor+ ((Plane*)allObjects[winIndex])->lambertShader(myray,PL);

			    //finalColor = finalColor+  ((Plane*)allObjects[winIndex])->getColor();
		    }


		    //finalColor = finalColor+  ((Plane*)allObjects[winIndex])->getColor();
		    if(mywinIndex!=-1 && softShadowFlag==0)
			    finalColor = finalColor+  Color(0,0,0);
		    if(cDirect==0){

			    finalColor =  finalColor*shadowColor;
                      }

		    if(!planeTextureEnabled){
                      finalColor = finalColor*(1-ks) + colorFromReflectedObject*ks;
                      return finalColor;

                    }

		    double X,Y,S0=10,S1=10;
		    // Edit the below line.. set S0 to 500 for just one picture.
		    S0=100;
		    S1=S0;
		    Point3f planeNorm =  ((Plane*)allObjects[winIndex])->normalVector;
		    planeNorm.Normalize();
		    Point3f crossVec = planeNorm^Point3f(0,0,1);
		    if(crossVec.x == 0 && crossVec.y == 0 && crossVec.z == 0 )
		    {

			    X =( (interSectionPoint - ((Plane*)allObjects[winIndex])->origin)%Point3f(1,0,0))/S0; 
			    Y =( (interSectionPoint - ((Plane*)allObjects[winIndex])->origin)%Point3f(0,1,0))/S1; 
		    }
		    crossVec = planeNorm^Point3f(0,1,0);
		    if(crossVec.x == 0 && crossVec.y == 0 && crossVec.z == 0 )
		    {

			    X =( (interSectionPoint - ((Plane*)allObjects[winIndex])->origin)%Point3f(1,0,0))/S0; 
			    Y =( (interSectionPoint - ((Plane*)allObjects[winIndex])->origin)%Point3f(0,0,1))/S1; 
		    }
		    crossVec = planeNorm^Point3f(1,0,0);
		    if(crossVec.x == 0 && crossVec.y == 0 && crossVec.z == 0 )
		    {

			    X =( (interSectionPoint - ((Plane*)allObjects[winIndex])->origin)%Point3f(0,1,0))/S0; 
			    Y =( (interSectionPoint - ((Plane*)allObjects[winIndex])->origin)%Point3f(0,0,1))/S1; 
		    }
		    double u =X- (int)X,v = Y - (int)Y;
		    if(u<0)	u = u+1;
		    if(v<0)	v = v+1;

		    //if(( (X>0 && X<1) && (Y>0 && Y<1) ) )   {                                  u = X; v = Y;
		    crossVec = planeNorm^Point3f(0,-1,0);
		    //if(crossVec.x == 0 && crossVec.y == 0 && crossVec.z == 0 )	// if its the roof plane, shade with this texture bumpMap
		    if(((Plane*)allObjects[winIndex])->id == "roof")
		    {
			    u = u*bumpImageWidth,v=v*bumpImageHeight;
			    int pixmapIndex = abs((int)v * bumpImageWidth + (int)u) * 3;
			    Color tmpColor;
			    tmpColor.red = (float)(bumpMap[pixmapIndex])/bumpmaxcolor;
			    tmpColor.green =(float)(bumpMap[pixmapIndex + 1])/bumpmaxcolor;
			    tmpColor.blue = (float)(bumpMap[pixmapIndex + 2])/bumpmaxcolor;
			   // cout<<"tmp roof ";tmpColor.printColor();cout<<endl;

			    finalColor = finalColor*tmpColor;

		//	cout<<"debug roof\n";
		    }
		    // crossVec = planeNorm^Point3f(0,0,-1);
		    //if(crossVec.x == 0 && crossVec.y == 0 && crossVec.z == 0 )	// if its the front plane, shade with this texture 
		     if(((Plane*)allObjects[winIndex])->id == "front")
		    {
			    u = u*projectionImageWidth1,v=v*projectionImageHeight1;
			    int pixmapIndex = abs((int)v * projectionImageWidth1 + (int)u) * 3;
			    Color tmpColor;
			    tmpColor.red = (float)(pixmap1[pixmapIndex])/maxcolor1;
			    tmpColor.green =(float)(pixmap1[pixmapIndex + 1])/maxcolor1;
			    tmpColor.blue = (float)(pixmap1[pixmapIndex + 2])/maxcolor1;
			   // cout<<"tmp front: ";tmpColor.printColor();cout<<endl;
			    finalColor = finalColor*tmpColor;
		//	cout<<"debug z=-1\n";
		    }
		    //crossVec = planeNorm^Point3f(0,1,0);
		    //if(crossVec.x == 0 && crossVec.y == 0 && crossVec.z == 0 )	// if its the floor plane, shade with this texture 
		     if(((Plane*)allObjects[winIndex])->id == "floor")
		    {
			    u = u*projectionImageWidth2,v=v*projectionImageHeight2;
			    int pixmapIndex = abs((int)v * projectionImageWidth2 + (int)u) * 3;
			    Color tmpColor;
			    tmpColor.red = (float)(pixmap2[pixmapIndex])/maxcolor2;
			    tmpColor.green =(float)(pixmap2[pixmapIndex + 1])/maxcolor2;
			    tmpColor.blue = (float)(pixmap2[pixmapIndex + 2])/maxcolor2;
			  //  cout<<"tmp floor: ";tmpColor.printColor();cout<<endl;
			    finalColor = finalColor*tmpColor;	
			 //   cout<<"debug floor\n";

		    }
		    //crossVec = planeNorm^Point3f(0,0,1);
		    //if(crossVec.x == 0 && crossVec.y == 0 && crossVec.z == 0 )	// if its the z=1 plane, shade with this texture 
		     if(((Plane*)allObjects[winIndex])->id == "back")
		    {
			    u = u*projectionImageWidth1,v=v*projectionImageHeight1;
			    int pixmapIndex = abs((int)v * projectionImageWidth1 + (int)u) * 3;
			    Color tmpColor;
			    tmpColor.red = (float)(pixmap1[pixmapIndex])/maxcolor1;
			    tmpColor.green =(float)(pixmap1[pixmapIndex + 1])/maxcolor1;
			    tmpColor.blue = (float)(pixmap1[pixmapIndex + 2])/maxcolor1;
			   // cout<<"tmp back: ";tmpColor.printColor();cout<<endl;
			    finalColor = finalColor*tmpColor;

		//	    cout<<"debug z=1\n";
		    }
		    
		     if(((Plane*)allObjects[winIndex])->id == "left")
		    {
			    u = u*projectionImageWidth1,v=v*projectionImageHeight1;
			    int pixmapIndex = abs((int)v * projectionImageWidth1 + (int)u) * 3;
			    Color tmpColor;
			    tmpColor.red = (float)(pixmap1[pixmapIndex])/maxcolor1;
			    tmpColor.green =(float)(pixmap1[pixmapIndex + 1])/maxcolor1;
			    tmpColor.blue = (float)(pixmap1[pixmapIndex + 2])/maxcolor1;
			   // cout<<"tmp back: ";tmpColor.printColor();cout<<endl;
			    finalColor = finalColor*tmpColor;

		//	    cout<<"debug z=1\n";
		    }
		    
		    if(((Plane*)allObjects[winIndex])->id == "right")
		    {
			    u = u*projectionImageWidth1,v=v*projectionImageHeight1;
			    int pixmapIndex = abs((int)v * projectionImageWidth1 + (int)u) * 3;
			    Color tmpColor;
			    tmpColor.red = (float)(pixmap1[pixmapIndex])/maxcolor1;
			    tmpColor.green =(float)(pixmap1[pixmapIndex + 1])/maxcolor1;
			    tmpColor.blue = (float)(pixmap1[pixmapIndex + 2])/maxcolor1;
			   // cout<<"tmp back: ";tmpColor.printColor();cout<<endl;
			    finalColor = finalColor*tmpColor;

		//	    cout<<"debug z=1\n";
		    }
	    }

	    if(allObjects[winIndex]->objectName=="genericObject")	
	    {
		    //cout<<"I do exist\n";
		    int numPlanes = ((GenericObject*)(allObjects[winIndex]))->triangles.faces.size(), faceindex =-1;
		    vector<Point3f> tmpVector;
		    //cout<<"numplanes: "<<numPlanes<<endl;
		    Point3f Ph = interSectionPoint;
		    //vector<FacePoint> vectPoints;



		    for(int i = 0; i < numPlanes; i++ )
		    {
			    ObjMeshFace triangleFace = ((GenericObject*)allObjects[winIndex])->triangles.faces[i];
			    Point3f triNorm(triangleFace.vertices[0].normal); 
			    triNorm.Normalize();
			    Plane triPlane( triNorm,triangleFace.vertices[0].pos,"triplane" );
			    //print(triPlane.origin);
			    Point3f P0 = triangleFace.vertices[0].pos;
			    Point3f P1 = triangleFace.vertices[1].pos;
			    Point3f P2 = triangleFace.vertices[2].pos;

			    //Point3f interPoint = triPlane.getIntersectionPoint(ray) , Ph = interPoint;

			    Point3f A0 = areaOfTriangle(P1, P2 , Ph);
			    Point3f A1 = areaOfTriangle(P2, P0 , Ph);
			    Point3f A2 = areaOfTriangle(P0, P1 , Ph);
			    Point3f A  = areaOfTriangle(P0, P1 , P2);

			    //cout<<"wx "<<A0.x/A.x<<endl;
			    //cout<<"wy "<<A0.y/A.y<<endl;
			    //cout<<"wz "<<A0.z/A.z<<endl;
			    double w0 = ( A0.x + A0.y + A0.z )/(A.x + A.y + A.z);
			    double w1 = ( A1.x + A1.y + A1.z )/(A.x + A.y + A.z);
			    double w2 = ( A2.x + A2.y + A2.z )/(A.x + A.y + A.z);

			    //printf("three ws: %f, %f, %f \n",w0,w1,w2);

			    // if in same plane
			    //if( w0 + w1 + w2 >= 0.99  && w0+w1+w2<=1.01 )
			    //{
			    // if inside the triangle
			    //cout<<"in the same plane!!\n";
			    double sumArea = A1.Length() + A2.Length() + A0.Length();

			    //if( (w0 >=0 && w0 <=1) && (w1 >=0 && w1 <=1) &&  (w2 >=0 && w2 <=1) )
			    if(sumArea >= 0.97*A.Length() && sumArea <= 1.03*A.Length())
			    {
				    //cout<<"inside the triangle!!\n";

				    faceindex = i;
				    //cout<<"face index of this point: ";print(interSectionPoint);cout<<endl;
				    Point3f centroid((triangleFace.vertices[0].pos + triangleFace.vertices[1].pos + triangleFace.vertices[2].pos)/3.0);
				  				    break;
			    }

			    //}
		    }
		    //cout<<endl;

		    //if(vectPoints.size() == 1)	faceindex = vectPoints[0].faceindex;
		    //if(vectPoints.size() == 2)
		    //	faceindex =( (Pe-vectPoints[0].point).Length() > (Pe-vectPoints[1].point).Length()) ? vectPoints[1].faceindex : vectPoints[0].faceindex;

		    finalColor = ((GenericObject*)allObjects[winIndex])->getColor();
		    if(!genericTextureEnabled){
			    //finalColor = finalColor*((GenericObject*)allObjects[winIndex])->goochShader(myray, PL);
			    finalColor= finalColor*(1-ks) + colorFromReflectedObject*ks;

                            //if(depth < numRecursion && eta>0.05)
                            //  finalColor = colorFromRefractedObject;
                            if(ks!=0 || eta!=0)
                            finalColor = (colorFromRefractedObject*(1-ks) + colorFromReflectedObject*ks);

                            return finalColor;

		    }
		    //if(( (X>0 && X<1) && (Y>0 && Y<1) ) )   {                                  u = X; v = Y;
		    /*
		       double u = ((GenericObject*)allObjects[winIndex])->triangles.faces[faceindex].vertices[0].texcoord.x;
		       double v = ((GenericObject*)allObjects[winIndex])->triangles.faces[faceindex].vertices[0].texcoord.y;

		       u = u*bumpImageWidth,v=v*bumpImageHeight;
		       int pixmapIndex = abs((int)v * bumpImageWidth + (int)u) * 3;
		       Color tmpColor;


		       tmpColor.red = (float)(bumpMap[pixmapIndex])/bumpmaxcolor;
		       tmpColor.green =(float)(bumpMap[pixmapIndex + 1])/bumpmaxcolor;
		       tmpColor.blue = (float)(bumpMap[pixmapIndex + 2])/bumpmaxcolor;

		       finalColor = finalColor*tmpColor;

		    //
		    if(faceindex ==0 )	finalColor = Color(1,1,1);
		    if(faceindex ==1 )	finalColor = Color(1,0,1);
		    if(faceindex ==2 )	finalColor = Color(1,1,0);
		    if(faceindex ==3 )	finalColor = Color(0,1,1);
		    */

		    //finalColor = finalColor*((GenericObject*)allObjects[winIndex])->goochShader(myray, PL);
		    //finalColor = finalColor* shadowColor;


		    for(int i = 0; i < numPlanes; i++ )
		    {
			    ObjMeshFace triangleFace = ((GenericObject*)allObjects[winIndex])->triangles.faces[i];
			    Point3f triNorm(triangleFace.vertices[0].normal); 
			    triNorm.Normalize();
			    Plane triPlane( triNorm,triangleFace.vertices[0].pos,"triplane" );
			    //print(triPlane.origin);
			    Point3f P0 = triangleFace.vertices[0].pos;
			    Point3f P1 = triangleFace.vertices[1].pos;
			    Point3f P2 = triangleFace.vertices[2].pos;

			    //Point3f interPoint = triPlane.getIntersectionPoint(ray) , Ph = interPoint;

			    Point3f A0 = areaOfTriangle(P1, P2 , Ph);
			    Point3f A1 = areaOfTriangle(P2, P0 , Ph);
			    Point3f A2 = areaOfTriangle(P0, P1 , Ph);
			    Point3f A  = areaOfTriangle(P0, P1 , P2);

			    //cout<<"wx "<<A0.x/A.x<<endl;
			    //cout<<"wy "<<A0.y/A.y<<endl;
			    //cout<<"wz "<<A0.z/A.z<<endl;
			    double w0 = ( A0.x + A0.y + A0.z )/(A.x + A.y + A.z);
			    double w1 = ( A1.x + A1.y + A1.z )/(A.x + A.y + A.z);
			    double w2 = ( A2.x + A2.y + A2.z )/(A.x + A.y + A.z);

			    //printf("three ws: %f, %f, %f \n",w0,w1,w2);

			    // if in same plane
			    //if( w0 + w1 + w2 >= 0.99  && w0+w1+w2<=1.01 )
			    //{
			    // if inside the triangle
			    //cout<<"in the same plane!!\n";
			    double sumArea = A1.Length() + A2.Length() + A0.Length();

			    //if( (w0 >=0 && w0 <=1) && (w1 >=0 && w1 <=1) &&  (w2 >=0 && w2 <=1) )
			    if(sumArea >= 0.97*A.Length() && sumArea <= 1.03*A.Length())
			    {
				    //cout<<"inside the triangle!!\n";

				    faceindex = i; 
				    double X,Y,S0=10,S1=10,Z;
				    S0=500;
				    S1=S0;


				    interSectionPoint =triNorm;
				    //interSectionPoint = interSectionPoint/((Sphere*)allObjects[winIndex])->radius;

				    X = interSectionPoint%Point3f(1,0,0);
				    Y = interSectionPoint%Point3f(0,1,0);
				    Z = interSectionPoint%Point3f(0,0,1);

				    /*
				       double u =X- (int)X,v = Y - (int)Y, w= Z-(int)Z;
				       if(u<0)	u = u+1;
				       if(v<0)	v = v+1;
				       if(w<0)	w = w+1;
				       */
				    double u =X- (int)X,v = Y - (int)Y;
				    if(u<0)	u = u+1;
				    if(v<0)	v = v+1;

				    //if(( (X>0 && X<1) && (Y>0 && Y<1) ) )   {                                  u = X; v = Y;
				    double u0 = triangleFace.vertices[0].texcoord.x, u1 = triangleFace.vertices[1].texcoord.x, u2 = triangleFace.vertices[2].texcoord.x ;
				    double v0 = triangleFace.vertices[0].texcoord.y, v1 = triangleFace.vertices[1].texcoord.y, v2 = triangleFace.vertices[2].texcoord.y ;

				    u = w0*u0 + w1*u1 + w2*u2;
				    v = w0*v0 + w1*v1 + w2*v2;
				    if(u<0)	u = u+1;
				    if(v<0)	v = v+1;
				    if(u>0)	u = u-1;
				    if(v>0)	v = v-1;


				    //if(u>1 || v>1|| u<0 || v<0)  
				    //  printf("absurd value: u: %f, v: %f", u,v);

				    u = u*projectionImageWidth,v=v*projectionImageHeight;
				    int pixmapIndex = abs((int)v * projectionImageWidth + (int)u) * 3;

				    //printf("psi: %f , theta: %f , u: %f , v: %f \n",psi, theta, u , v);
				    //cout<<(int)( (Y * projectionImageWidth + X) * 3 )<<endl; 
				    //cout<< (float)pixmap[pixmapIndex]<<endl;
				    finalColor.red = (float)(pixmap[pixmapIndex])/maxcolor;
				    finalColor.green =(float)(pixmap[pixmapIndex + 1])/maxcolor;
				    finalColor.blue = (float)(pixmap[pixmapIndex + 2])/maxcolor;


				    //cout<<"face index of this point: ";print(interSectionPoint);cout<<endl;
				    break;
			    }

			    //}
		    }

		    finalColor = finalColor*((GenericObject*)allObjects[winIndex])->goochShader(myray, PL);
		    //finalColor = finalColor * shadowColor;


		    }
		     if(allObjects[winIndex]->objectName=="genericObject"){
		     	finalColor = finalColor * (finalColor*(1-ks) + colorFromReflectedObject*ks);
		     }
		     else {
		     	finalColor = (finalColor*(1-ks) + colorFromReflectedObject*ks);
		     }
		    
                    if(depth < numRecursion && eta>0.05)
                      finalColor =  colorFromRefractedObject;
    return finalColor;
}

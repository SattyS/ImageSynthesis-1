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
SpotLight spotLight(Point3f(0,70,45),Point3f(0,0,1),60.0/180.0);
//Enable each variable to enable textures on them
bool sphereTextureEnabled=false, genericTextureEnabled=false,planeTextureEnabled=false, textureRefractionMapEnabled=false, environmentFlag = true;
bool glossyEnabled=0;
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
int numRecursion =1;

Point3f geoCentre(0,0,0);

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
        initGeo();
	
	RGBType *pixels= new RGBType[n];
	int index=0,M=1,N=M; 
	int Sx=10,winIndex=0;
	int Sy=(Sx*Ymax)/Xmax;
	float x,y;
	
	Point3f Vview(0,-10,45),Vup(0,1,0);	// point the view vector to focus on a particular point from Pe
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

	Sphere sphere1(Point3f(0,-10,45),10, Color(1,1,1),1,0 , 0);
	Sphere sphere2(Point3f(15,0,110),17, Color(1,0.2,0.3),2,0.9, 0);
	Sphere sphere3(Point3f(0,0,85),16, Color(0.1,0.5,1),3,0.9, 0);
	//Sphere sphere4(Point3f(1,1,6),3, Color(0,0.5,1),4);
	
	Plane plane1(Point3f(0,-1,0), Point3f(0,40,0), Color(1,1,1), "roof", 0,0);
	Plane plane2(Point3f(-1,0,0), Point3f(60,0,0), Color(1,0,0), "left", 0,0);
	Plane plane3(Point3f(0,0,-1), Point3f(0,0,60), Color(1,1,1), "front",0,0);
	Plane plane4(Point3f(0,0,1), Point3f(0,0,-80), Color(0,0,0), "back", 0,0);
	Plane plane5(Point3f(0,1,0), Point3f(0,-20,0), Color(0.1,0.15,0.51), "floor",0,0);
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
	//allObjects.push_back(dynamic_cast<Object*>(&sphere2));
	//allObjects.push_back(dynamic_cast<Object*>(&sphere3));
	//allObjects.push_back(dynamic_cast<Object*>(&sphere4));
	
	//allObjects.push_back(dynamic_cast<Object*>(&plane1));
	//allObjects.push_back(dynamic_cast<Object*>(&plane2));
	//allObjects.push_back(dynamic_cast<Object*>(&plane3));
	//allObjects.push_back(dynamic_cast<Object*>(&plane4));
	allObjects.push_back(dynamic_cast<Object*>(&plane5));
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
        int softShadowFlag=1;
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
                        
                        // comment the following line for ENVIRONMENT MAP:  
                        Point3f refLectedRayDirection =myray.direction ;refLectedRayDirection.Normalize();
                        Color colorFromReflectedObject = finalColor;

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

                          finalColor =  colorFromReflectedObject;



                        }



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
    geoCentre = interSectionPoint;
    //for(int i=0;i<geodesicDome.size();i++)
      //geodesicDome[i] = geodesicDome[i] + interSectionPoint;

    int flag=0;
    Point3f NH ;
    finalColor = Color(0,0,0);

    Color colorFromRefractedObject(0,0,0);
    double ks=allObjects[winIndex]->KS;
    double eta=allObjects[winIndex]->findEta(interSectionPoint);
    float rnd1 =-3 + 6*(rand()/float(RAND_MAX)),rnd2 = -1 + 2*(rand()/float(RAND_MAX)),rnd3 = -1 + 2*(rand()/float(RAND_MAX)) ;
    //rnd1 = pow(rnd1,4);
    //rnd2 = pow(rnd2,4);
    //rnd3 = pow(rnd3,4);
    //printf(" rnd1,2,3: %f, %f, %f\n ",rnd1,rnd2,rnd3);

     TotalfinalColor = Color(0,0,0);


     // write the geodesic sphere and shadow equations!


    finalColor = Color(0,0,0);

int numRaysOutside = geodesicDome.size(), BigConst=1000000;

     for(int qwe = 0 ; qwe < geodesicDome.size() ; qwe++ )
     {

       Point3f geoPoint = geodesicDome[qwe];

       Point3f geoRayDir = geodesicDome[qwe];
       //Point3f geoRayDir = (geoPoint - interSectionPoint);
       geoRayDir.Normalize();


       float costeta = geoRayDir%allObjects[winIndex]->getNormal(interSectionPoint);
       if(costeta >1.0 || costeta <0.0) // check if it should shoot a ray
       {
         //cout<<"geor: ";print(geoRayDir);
         //cout<<"normal:";print(allObjects[winIndex]->getNormal(interSectionPoint));
         //cout<<"costeta"<<costeta<<endl;

         numRaysOutside--;
         continue;
       }
        
       PL = interSectionPoint + BigConst*geoRayDir;

       if(allObjects[winIndex]->objectName=="Sphere" )
       {
       
       }






    //int numLights = lights.size();
    //print(PL);
    
     
    Ray jujuRay(rayStart , PL-rayStart);

    flag=0;
    tmp.clear();
    
    
    
    
    // change PL here to environment map
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
            else if(allObjects[j]->objectName=="Plane" )
	    {
		    shadowDist +=abs( ( ((Plane*)allObjects[j])->normalVector % PL) -((Plane*)allObjects[j])->origin.y )  ;
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
		    //softShadowFlag = 1;
                    ratio = 0;


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
			    //finalColor = finalColor + shadowColor*(((Sphere*)allObjects[winIndex])->phongShader(myray,PL));
			    finalColor = finalColor + shadowColor*(((Sphere*)allObjects[winIndex])->getColor());
                            //cout<<"eta: "<<eta<<endl;
                            //cout<<"colorFromRefractedObject: "; colorFromRefractedObject.printColor();cout<<endl;
			    //if(ks!=0 || eta!=0)
                            //finalColor =(colorFromRefractedObject*(1-ks) + colorFromReflectedObject*ks);
                            //cout<<"finalColor: "; finalColor.printColor();cout<<endl;
			    /*
                             * float myalpha = 0.3;
                            if(depth < numRecursion && eta>0.05){
                              finalColor = finalColor*myalpha +colorFromRefractedObject*(1-myalpha);
                              finalColor = colorFromRefractedObject;
                            }
                            */

			    //return finalColor;
		    }
                    else{
		    interSectionPoint =( interSectionPoint - ((Sphere*)allObjects[winIndex])->center);
		    interSectionPoint = interSectionPoint/((Sphere*)allObjects[winIndex])->radius;

		    X = interSectionPoint%Point3f(1,0,0);
		    Y = interSectionPoint%Point3f(0,1,0);
		    Z = interSectionPoint%Point3f(0,0,1);
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
	    }

            }
	    if(mywinIndex!=-1 && softShadowFlag==0)
		    finalColor = finalColor+  Color(0,0,0);

            if(allObjects[winIndex]->objectName=="Plane")	
            {
		    if(spotlightEnabled)
			    finalColor = finalColor+ shadowColor*(((Plane*)allObjects[winIndex])->phongShader(myray,PL))*alpha;
		    else
		    {
			    //finalColor = finalColor+ shadowColor*((Plane*)allObjects[winIndex])->lambertShader(myray,PL);
			    finalColor = finalColor+  shadowColor*((Plane*)allObjects[winIndex])->getColor();
		    }

		    //finalColor = finalColor+  ((Plane*)allObjects[winIndex])->getColor();
		    if(mywinIndex!=-1 && softShadowFlag==0)
			    finalColor = finalColor+  Color(0,0,0);
	
		    if(!planeTextureEnabled){
                      //finalColor = finalColor*(1-ks) + colorFromReflectedObject*ks;
                      //return finalColor;

                    }
                  else
                  {
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

	    }
              }   // end of geodesic light loop


                        //finalColor = finalColor/(geodesicDome.size()/2.0);
                        finalColor = finalColor/float(numRaysOutside);
                        //cout<<numRaysOutside<<"\n";
                
                  return finalColor;
            }

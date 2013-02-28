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

using namespace std;
using namespace cy;

int main (int argc, char const* argv[])
{
	Point3f cr(-4.152,-0.534,-14.12);
	cr.Normalize();
	//float dot= p%q;
	cout<<cr.x<<" "<<cr.y<<" "<<cr.z;cout<<endl;
	//cout<<dot;cout<<endl;
	float r=29.6,c=0;
	Point3f pe(-0.447,-3.725,1.465),pc(-6,-4,-12),npe(4.0/30.0,1.0/3.0,28.0/30.0);
	
	c=(pe-pc)%(pe-pc) - r*r;
	printf("c: %f\n",c);//cout<<c<<endl;
	float b=npe % (pc-pe); 
	printf("b: %f\n",b);cout<<b<<endl;
	printf("delta = b2 - c = %f\n",b*b-c);
	float delta=b*b-c, th=b-sqrt(delta);
	printf("th = %f\n",th);
	Point3f ph=pe+(npe*th), nh=(ph-pc)/r;
	
	printf("ph = (%f,%f,%f) \n",ph.x,ph.y,ph.z);
	
	ph=(pe-pc);ph.Normalize();
	printf("ph = (%f,%f,%f) \n",ph.x,ph.y,ph.z);
	
	printf("nh = (%f,%f,%f) \n",nh.x,nh.y,nh.z);
	//npe=(pe^pc).Normalize();
	
	ph=Point3f(16.0,3.0,3.0);ph.Normalize();
	printf("ph = (%f,%f,%f) \n",ph.x,ph.y,ph.z);
	Point3f pL(20,7,10),n(0,0,1);
	
	printf("cos0 = %f\n", (pL%n)/pL.Length());
	
	
	//Point3f R=npe;
	return 0;
}

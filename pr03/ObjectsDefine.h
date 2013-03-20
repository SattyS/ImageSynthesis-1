#include "saveAndRay.h"

using namespace cy;
using namespace std;
#define print(p) printf("(%f,%f,%f) \n",p.x,p.y,p.z);
#define pb(a) push_back(a)
extern Point3f Pe;
//extern Point3f PL;

/*************************************** OBJECT CLASS **********************************************************************/
class Object 
{
	public:
		string objectName;
		Object(){}
		virtual Point3f getIntersectionPoint(Ray ray){return Point3f(0,0,0);}
		virtual bool isEyeOutside(Point3f Pe){return false;}
		virtual Color getColor(){return Color(0,0,0);}
		virtual Color lambertShader(Ray myray,Point3f PL){return Color(0,0,0);}
		virtual Color goochShader(Ray ray,Point3f PL){return Color(0,0,0);}

		virtual Color phongShader(Ray ray,Point3f PL){return Color(0,0,0);}

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

			//cout<<"value t: "<<t<<" ";
			if (t<0)	return ret;
			ret = ray.origin + t*ray.direction;
			//print(ret);
			return ret;
		}	
		virtual Color lambertShader(Ray myray,Point3f PL)
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
			//cout<< "cos: "<<(c)<<endl;
			if((c)<0)	//c=0.0;
			{
				//finalColor = ambient_color;
				//return finalColor;
				c=0;
			}
			finalColor = diffused_color*c + ambient_color*(1-c);
			//s=pow(s,alpha);
			//cout<< "c: "<<c<<endl;
			// Gooch shading
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
		virtual Color goochShader(Ray myray,Point3f PL)
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

		virtual Color phongShader(Ray myray,Point3f PL)
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
			finalColor = lambertShader(myray,PL);

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
	double getTwoDelta(Ray ray)
	{
		//cout<<"assuming that eye is outside the sphere \n";
		double c=(center-ray.origin)%(center-ray.origin)-(radius*radius),b=(center-ray.origin)%ray.direction;
		double delta=b*b-c;
		if(delta<0)	{
			//cout<<"delta<0\n";
			return 0.0;
		}
		return 2*sqrt(delta);
	
	}
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
	virtual Color lambertShader(Ray myray,Point3f PL)
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
	virtual Color goochShader(Ray myray,Point3f PL)
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

	virtual Color phongShader(Ray myray,Point3f PL)
	{
		//========================================== L & G =======================================
		Color spColor ;
		Point3f Ph,Pc, nh,nlh;
		double cosTheta=0,c;
		double alpha,s,s0,delta;
		Point3f v,r;
		Color ambient_color(0,0,0),diffused_color,finalColor,specular_color(1,1,1);

		finalColor = lambertShader(myray,PL );
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



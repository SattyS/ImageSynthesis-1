#include "saveAndRay.h"
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdio.h>

#define TOKEN_VERTEX_POS "v"
#define TOKEN_VERTEX_NOR "vn"
#define TOKEN_VERTEX_TEX "vt"
#define TOKEN_FACE "f"


using namespace cy;
using namespace std;
#define print(p) printf("(%f,%f,%f) \n",p.x,p.y,p.z);
#define pb(a) push_back(a)
extern Point3f Pe;
//extern Point3f PL;

/**
 *  * The MIT License
 *   * 
 *    * Copyright (c) 2010 Wouter Lindenhof (http://limegarden.net)
 *     * 
 *      * Permission is hereby granted, free of charge, to any person obtaining a copy
 *       * of this software and associated documentation files (the "Software"), to deal
 *        * in the Software without restriction, including without limitation the rights
 *         * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *          * copies of the Software, and to permit persons to whom the Software is
 *           * furnished to do so, subject to the following conditions:
 *            * 
 *             * The above copyright notice and this permission notice shall be included in
 *              * all copies or substantial portions of the Software.
 *               * 
 *                * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *                 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *                  * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *                   * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *                    * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *                     * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *                      * THE SOFTWARE.
 *                       */
/*************************************** OBJECT CLASS **********************************************************************/
class Object 
{
	public:
		string objectName;
		float KS;
		Object(){KS=0;}
		virtual Point3f getIntersectionPoint(Ray ray){return Point3f(0,0,0);}
		virtual bool isEyeOutside(Point3f Pe){return false;}
		virtual Color getColor(){return Color(0,0,0);}
		virtual Point3f getNormal(Point3f P){}
		virtual Color lambertShader(Ray myray,Point3f PL){return getColor();}
		virtual Color goochShader(Ray ray,Point3f PL){return getColor();}

		virtual Color phongShader(Ray ray,Point3f PL){return getColor();}

};

Point3f areaOfTriangle(Point3f P0, Point3f P1, Point3f P2)
{
	// returns an area vector
	Point3f v1 = P1 - P0 , v2 = P2 - P0;
	Point3f v0 = v1^v2;
	//double area = v0.Length()/2.0;
	return v0/2.0;
}

struct ObjMeshVertex{
	Point3f pos;
	Point2f texcoord;
	Point3f normal;
};

/* This is a triangle, that we can render */
struct ObjMeshFace{
	ObjMeshVertex vertices[3];
};

/* This contains a list of triangles */
struct ObjMesh{
	std::vector<ObjMeshFace> faces;
};

/* Internal structure */
struct _ObjMeshFaceIndex{
	int pos_index[3];
	int tex_index[3];
	int nor_index[3];
};

 

/*************************************** OBJECT CLASS **********************************************************************/
class Plane: public Object
{
	public:
		
		Point3f normalVector;
		Point3f origin;
		Color color;
		Plane(){normalVector = Point3f(0,0,1);origin=Point3f(0,0,0);color=Color(1,0.5,0.2);objectName="Plane";KS=0;}
		Plane(Point3f r, Point3f orig){normalVector = r; origin=orig;	objectName="Plane";normalVector.Normalize();KS=0;}
		Plane(Point3f r, Point3f orig,Color c){normalVector = r; origin=orig; color = c; objectName="Plane";normalVector.Normalize();KS=0;}
		virtual bool isEyeOutside(Point3f Pe)
		{
			double c=(Pe-origin)%(normalVector);
			
			if(c>0)		//Eye is outside the sphere 
				return true;
			else 	return false;
		}
		virtual Color getColor(){	return color;}
		virtual Point3f getNormal(Point3f P){	return normalVector;}
		virtual Point3f getIntersectionPoint(Ray ray)
		{
			Point3f ret(-1,-1,-1);
			if( !isEyeOutside(ray.origin )  )	return ret;
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


/* Call this function to load a model, only loads triangulated meshes */
ObjMesh LoadObjMesh(std::string filename)
{
	ObjMesh myMesh;
	std::vector<Point3f>           positions;
	std::vector<Point2f>           texcoords;
	std::vector<Point3f>           normals;
	std::vector<_ObjMeshFaceIndex>  faces;
	/**
	  Load file, parse it
	 *           * Lines beginning with: 
	 *                * '#'  are comments can be ignored
	 *                     * 'v'  are vertices positions (3 floats that can be positive or negative)
	 *                          * 'vt' are vertices texcoords (2 floats that can be positive or negative)
	 *                               * 'vn' are vertices normals   (3 floats that can be positive or negative)
	 *                                    * 'f'  are faces, 3 values that contain 3 values which are separated by / and <space>
	 *                                         */

	std::ifstream filestream;
	filestream.open(filename.c_str());

	std::string line_stream;	// No longer depending on char arrays thanks to: Dale Weiler
	while(std::getline(filestream, line_stream)){	
		std::stringstream str_stream(line_stream);
		std::string type_str;
		str_stream >> type_str;
		if(type_str == TOKEN_VERTEX_POS){
			Point3f pos;
			str_stream >> pos.x >> pos.y >> pos.z;
			positions.push_back(pos);
		}else if(type_str == TOKEN_VERTEX_TEX){
			Point2f tex;
			str_stream >> tex.x >> tex.y;
			texcoords.push_back(tex);
		}else if(type_str == TOKEN_VERTEX_NOR){
			Point3f nor;
			str_stream >> nor.x >> nor.y >> nor.z;
			normals.push_back(nor);
		}else if(type_str == TOKEN_FACE){
			_ObjMeshFaceIndex face_index;
			char interupt;
			for(int i = 0; i < 3; ++i){
				str_stream >> face_index.pos_index[i] >> interupt 
					>> face_index.tex_index[i]  >> interupt 
					>> face_index.nor_index[i];
			}
			faces.push_back(face_index);
		}
	}
	// Explicit closing of the file 
	filestream.close();
//cout<<"just before :for: in LoadObjMesh function... \n";
	for(size_t i = 0; i < faces.size(); ++i){
		ObjMeshFace face;
		//cout<<"in: for1\n";
		for(size_t j = 0; j < 3; ++j){
		//cout<<"\t \t in: for2\n";
			face.vertices[j].pos        = positions[faces[i].pos_index[j] - 1];
			//cout<<"\t \t out: for2\n";
			face.vertices[j].texcoord   = texcoords[faces[i].tex_index[j] - 1];
			//cout<<"\t \t out: for2\n";
			face.vertices[j].normal     = normals[faces[i].nor_index[j] - 1];
			//cout<<"\t \t out: for2\n";
		}
		//cout<<"out: for1\n";
		myMesh.faces.push_back(face);
	}

cout<<"exiting the LoadObjMesh function... \n";
	return myMesh;
}

/*************************************** OBJECT for .obj files CLASS **********************************************************************/

class GenericObject:public Object
{
	public:
		ObjMesh triangles;
		Color color;
		//objLoader *objData; 


		GenericObject()
		{
			cout<<"constructor str filename\n";
			//objData = new objLoader();
			//objData->load("tetrahedron.obj");
			objectName = "genericObject";
			triangles = LoadObjMesh("tetrahedron.obj");
			color = Color(1,1,1);KS=0;
		}
		GenericObject(ObjMesh tri,float ks)
		{
			cout<<"constructor str filename\n";
			triangles = tri;		
			color = Color(1,0,0);
			objectName = "genericObject";KS=ks;
		}
		GenericObject(string fn, float ks)
		{
			cout<<"Entering LoadObjMesh function... \n";
			//objData = new objLoader();
			//objData->load(fn);
			triangles = LoadObjMesh(fn);
			cout<<"constructor str filename\n";
			objectName = "genericObject";
			color = Color(1,1,1);KS = ks;
		}
		virtual bool isEyeOutside(Point3f Pe)
		{
			int numPlanes = triangles.faces.size();
			bool val=true, flag =true;
			for(int i = 0; i < numPlanes; i++ )
			{
			
				ObjMeshFace triangleFace = triangles.faces[i];
				//obj_face *o = objData->faceList[i];

				Point3f triPlaneOrig, triPlaneNormal;
				//triPlaneOrig.x = objData->vertexList[ o->vertex_index[0] ];
				//triPlaneOrig.y = objData->vertexList[ o->vertex_index[1] ];
				//triPlaneOrig.z = objData->vertexList[ o->vertex_index[2] ];

				
				Plane triPlane( triangleFace.vertices[0].normal  , triangleFace.vertices[0].pos );
				if(i==0)	{ val = triPlane.isEyeOutside(Pe); continue; }

				if(val == triPlane.isEyeOutside(Pe)){
					flag &= true;
				}
				else{
					flag &= false;
					val = triPlane.isEyeOutside(Pe);
				}
				//cout<<"value of val,flag"<<val<<flag<<endl;
			}

			if(flag)	return false;

			return true;
		}

		virtual Point3f getNormal(Point3f Ph)
		{
			int numPlanes = triangles.faces.size(), faceindex =-1;
			vector<Point3f> tmpVector;Point3f nh;
			    //cout<<"numplanes: "<<numPlanes<<endl;
			    //vector<FacePoint> vectPoints;

			for(int i = 0; i < numPlanes; i++ )
			{
				ObjMeshFace triangleFace = triangles.faces[i];
				Point3f triNorm(triangleFace.vertices[0].normal); 
				triNorm.Normalize();
				Plane triPlane( triNorm,triangleFace.vertices[0].pos );
				//print(triPlane.origin);
				Point3f P0 = triangleFace.vertices[0].pos;
				Point3f P1 = triangleFace.vertices[1].pos;
				Point3f P2 = triangleFace.vertices[2].pos;
				//Point3f interPoint = triPlane.getIntersectionPoint(ray) , Ph = interPoint;
				Point3f A0 = areaOfTriangle(P1, P2 , Ph);
				Point3f A1 = areaOfTriangle(P2, P0 , Ph);
				Point3f A2 = areaOfTriangle(P0, P1 , Ph);
				Point3f A  = areaOfTriangle(P0, P1 , P2);
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
						break;
					}
				//}
			}
			//cout<<endl;
			// ========================== Lambert and Gooch =========================================	
			//Pc=origin;
			nh=triangles.faces[faceindex].vertices[0].normal;
			return nh;
		}
		virtual Color getColor()
		{
			return color;}

		virtual Point3f getIntersectionPoint(Ray ray)
		{
			int numPlanes = triangles.faces.size(), index =-1;
			vector<Point3f> tmpVector;
            //cout<<"numplanes: "<<numPlanes<<endl;

			if( !isEyeOutside(ray.origin )  )	return Point3f(-1,-1,-1);
			for(int i = 0; i < numPlanes; i++ )
			{
				ObjMeshFace triangleFace = triangles.faces[i];
				Point3f triNorm(triangleFace.vertices[0].normal); 
				triNorm.Normalize();
				Plane triPlane( triNorm,triangleFace.vertices[0].pos );
				//print(triPlane.origin);
				Point3f P0 = triangleFace.vertices[0].pos;
				Point3f P1 = triangleFace.vertices[1].pos;
				Point3f P2 = triangleFace.vertices[2].pos;

				Point3f interPoint = triPlane.getIntersectionPoint(ray) , Ph = interPoint;

				Point3f A0 = areaOfTriangle(P1, P2 , Ph);
				Point3f A1 = areaOfTriangle(P2, P0 , Ph);
				Point3f A2 = areaOfTriangle(P0, P1 , Ph);
				Point3f A  = areaOfTriangle(P0, P1 , P2);

				double w0 = ( A0.x + A0.y + A0.z )/(A.x + A.y + A.z);
				double w1 = ( A1.x + A1.y + A1.z )/(A.x + A.y + A.z);
				double w2 = ( A2.x + A2.y + A2.z )/(A.x + A.y + A.z);

				//printf("three ws: %f, %f, %f \n",w0,w1,w2);

				// if in same plane
				if( w0 + w1 + w2 >= 0.95 && w0 + w1 + w2 <=1.05 )
				{
					// if inside the triangle
					//cout<<"in the same plane!!\n";
					if( (w0 >=0 && w0 <=1) && (w1 >=0 && w1 <=1) &&  (w2 >=0 && w2 <=1) )
					{
					    //cout<<"inside the triangle!!\n";
						tmpVector.push_back(interPoint);
					}
				 }
			}
			double minDist = INT_MAX, len=0.0;
			//cout<<tmpVector.size()<<endl;
			for(int i = 0; i < tmpVector.size(); i++)
			{
				// if plane doesn't have any intersection point
				//if(tmpVector[i].x == -1 && tmpVector[i].y == -1 && tmpVector[i].z == -1)	continue;

				len=(tmpVector[i]-ray.origin).Length();
				if(len < minDist)
				{
					minDist=len;
					index=i;
				}
			}
			if (index == -1)	return Point3f(-1,-1,-1);
			//cout<<"generic interpoint: i was called! \n";
			//print(tmpVector[index]);
			return tmpVector[index];

		}

                double twoDist(Ray ray)
                {
                  int numPlanes = triangles.faces.size(), index =-1;
                  vector<Point3f> tmpVector;
            //cout<<"numplanes: "<<numPlanes<<endl;

			for(int i = 0; i < numPlanes; i++ )
			{
				ObjMeshFace triangleFace = triangles.faces[i];
				Point3f triNorm(triangleFace.vertices[0].normal); 
				triNorm.Normalize();
				Plane triPlane( triNorm,triangleFace.vertices[0].pos );
				//print(triPlane.origin);
				Point3f P0 = triangleFace.vertices[0].pos;
				Point3f P1 = triangleFace.vertices[1].pos;
				Point3f P2 = triangleFace.vertices[2].pos;

				Point3f interPoint = triPlane.getIntersectionPoint(ray) , Ph = interPoint;

				Point3f A0 = areaOfTriangle(P1, P2 , Ph);
				Point3f A1 = areaOfTriangle(P2, P0 , Ph);
				Point3f A2 = areaOfTriangle(P0, P1 , Ph);
				Point3f A  = areaOfTriangle(P0, P1 , P2);

				double w0 = ( A0.x + A0.y + A0.z )/(A.x + A.y + A.z);
				double w1 = ( A1.x + A1.y + A1.z )/(A.x + A.y + A.z);
				double w2 = ( A2.x + A2.y + A2.z )/(A.x + A.y + A.z);

				//printf("three ws: %f, %f, %f \n",w0,w1,w2);

				// if in same plane
				if( w0 + w1 + w2 >= 0.95 && w0 + w1 + w2 <=1.05 )
				{
					// if inside the triangle
					//cout<<"in the same plane!!\n";
					if( (w0 >=0 && w0 <=1) && (w1 >=0 && w1 <=1) &&  (w2 >=0 && w2 <=1) )
					{
					    //cout<<"inside the triangle!!\n";
						tmpVector.push_back(interPoint);
					}
				 }
			}
			double minDist = INT_MAX, len=0.0;
			//cout<<tmpVector.size()<<endl;

                        if(tmpVector.size() ==2 )
                        {
                          return (tmpVector[0] - tmpVector[1]).Length();  
                        
                        }
                        //else
                        //  cout<<"something is wrong ... check!\n";
			return 0.0;



                }
		virtual Color lambertShader(Ray myray,Point3f PL){
		
			return getColor();}
		virtual Color goochShader(Ray ray,Point3f PL){
			
			Color spColor ;
			Point3f Ph,Pc, nh,nlh;
			double cosTheta=0,c;
			double alpha,s,s0,delta;
			Point3f v,r;
			Color ambient_color(0,0,0),diffused_color,finalColor,specular_color(1,1,1);
			//for(int i=0; i < allObjects.size(); i++)
			//{
			spColor = color;
			Ph=getIntersectionPoint(ray);
			//cout<<"I do exist\n";
			int numPlanes = triangles.faces.size(), faceindex =-1;
			    vector<Point3f> tmpVector;
			    //cout<<"numplanes: "<<numPlanes<<endl;
			    //vector<FacePoint> vectPoints;

			    

			for(int i = 0; i < numPlanes; i++ )
			{
				ObjMeshFace triangleFace = triangles.faces[i];
				Point3f triNorm(triangleFace.vertices[0].normal); 
				triNorm.Normalize();
				Plane triPlane( triNorm,triangleFace.vertices[0].pos );
				//print(triPlane.origin);
				Point3f P0 = triangleFace.vertices[0].pos;
				Point3f P1 = triangleFace.vertices[1].pos;
				Point3f P2 = triangleFace.vertices[2].pos;

				//Point3f interPoint = triPlane.getIntersectionPoint(ray) , Ph = interPoint;

				Point3f A0 = areaOfTriangle(P1, P2 , Ph);
				Point3f A1 = areaOfTriangle(P2, P0 , Ph);
				Point3f A2 = areaOfTriangle(P0, P1 , Ph);
				Point3f A  = areaOfTriangle(P0, P1 , P2);
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
						break;
					}
					
				//}
			}
			//cout<<endl;


			// ========================== Lambert and Gooch =========================================	
			diffused_color = color;
			//Pc=origin;

			nh=triangles.faces[faceindex].vertices[0].normal;

			nlh=(PL - Ph);nlh.Normalize();
			cosTheta = nh % nlh;

			// try different values with c and cosTheta
			c = cosTheta;
			c = (cosTheta+1.0)/2.0;
			if(c<0)	c=0.0;
			finalColor = diffused_color*c + ambient_color*(1-c);
			return finalColor;
			
			//return getColor();
		}

		virtual Color phongShader(Ray ray,Point3f PL)
		{
		
		
		return getColor();}
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
		objectName = "Sphere";KS=0;
	}
	Sphere(Point3f c, double rad,Color col,int i, float ks)
	{
		radius= rad;
		id=i;
		center=c;
		color=col;
		objectName = "Sphere";KS = ks;
	}
	Sphere(Point3f c, double rad, int i,float ks)
	{
		id=i;
		radius= rad;
		center=c;
		color=Color(1,1,1);
		objectName = "Sphere";KS = ks;
	}
	
	virtual Point3f getNormal(Point3f point)
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
		Point3f ret(-1,-1,-1);
		//cout<<"assuming that eye is outside the sphere \n";
		
		if( !isEyeOutside(ray.origin )  )	return ret;
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



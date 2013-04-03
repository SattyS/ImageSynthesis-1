/////////////////////////////////////////////////////////////////////////
//////////////////////////// Light Classes //////////////////////////////
/////////////////////////////////////////////////////////////////////////

class Light
{
public:
  Color Lcolor;
  float intensity;
  int type;
  // type=1: pointlight
  // type=2: spotlight
  // type=3: directional light
  // type=4: arealight


  V3 position;
  V3 direction;

  float cosAngle; // light coneAngle, onle apply for spotligt

  float Lsx, Lsy; // light size, only apply for arealight
  int LM, LN; //light sample size, only apply for arealight
  V3 n0; //local coordinate system, onle apply for arealight



  Light(int _t, Color c, float _in, V3 pt, V3 d, float _ca, float _lsx, float _lsy, int _lm, int _ln, V3 v){
    // lightType, lightColor, intensity, position, direction, coneAngle, lightSize, lightSample, lightLocalCoor
    type = _t; 
    Lcolor = c;
    intensity = _in;

    position = pt;
    direction = d;
    cosAngle = _ca;

    Lsx = _lsx;
    Lsy = _lsy;
    LM = _lm;
    LN = _ln;
    n0 = v;

  }


  V3 getn1(){ 
     V3 n1(direction^n0); 
     n1.Normalize();
     return n1;
  };

  virtual  V3 getShadowRay(V3 Pt){};
  virtual  float computShadow(V3 Pt, int index, V3 Ray, int rand_x, int rand_y){};
  virtual  Color addShadowColor(Color shadColor, float sh){};

};




//////////////////////////// POINTLIGHT /////////////////////////////////

class pointLight: public Light
{
public:

  //V3 position;

  //construtor
  pointLight (V3 pt, Color c, float _in):Light(1, c, _in, pt, pt, 1.0, 0, 0, 0, 0, pt)
  {
    //position = pt
    //lightcolor = c
    //lightIntensity = _in
  }



  V3 getShadowRay(V3 Pt){

    V3 shadowRay(position - Pt);
    Dlh = shadowRay.Length();
    shadowRay.Normalize();
    return shadowRay;

  };

  float computShadow(V3 Pt, int index, V3 Ray, int rand_x, int rand_y);
  Color addShadowColor(Color shadColor, float sh){};

private:

  float Dlh;

};



//////////////////////////// SPOTLIGHT /////////////////////////////////

class spotLight: public Light
{
public:

  //V3 position;
  //V3 direction;
  //float cosAngle;

  //construtor
  spotLight (V3 pt, V3 d, float _ca, Color c, float _in):Light(2, c, _in, pt, d, _ca, 0, 0, 0, 0, pt)
  {
    //position = pt
    //direction = d
    //cosAngle = ca 
    //lightcolor = c
    //lightIntensity = _in
  }

  V3 getShadowRay(V3 Pt){

    V3 shadowRay(position - Pt);
    Dlh = shadowRay.Length();
    shadowRay.Normalize();
    return shadowRay;

  };

  float computShadow(V3 Pt, int index, V3 Ray, int rand_x, int rand_y);
  Color addShadowColor(Color shadColor, float sh){};


private:

  float Dlh;

};




//////////////////////////// DIRECTIONAL_LIGHT /////////////////////////////////

class directionalLight: public Light
{
public:

  //V3 direction;

  //construtor
  directionalLight (V3 d, Color c, float _in):Light(3, c, _in, d, d, 0, 0, 0, 0, 0, d)
  {
    //direction = d
    //lightcolor = c
    //lightIntensity = _in
  }


  V3 getShadowRay(V3 Pt){

    V3 shadowRay((-1)*direction);
    return shadowRay;

  };


  float computShadow(V3 Pt, int index, V3 Ray, int rand_x, int rand_y);
  Color addShadowColor(Color shadColor, float sh){};


private:

  float Dlh;

};




//////////////////////////// AREALIGHT /////////////////////////////////

class areaLight: public Light
{
public:

  //construtor
  areaLight (V3 pt, V3 nz, V3 nx, float _lsx, float _lsy, int _lm, int _ln, Color c, float _in):Light(4, c, _in, pt, nz, 1.0, _lsx, _lsy, _lm, _ln, nx)
  {
    //position = pt
    //normal = nz
    //UpVector = nx
    //LightSize = _lsx, _lsy
    //SampleNum = _lm, _ln
    //lightColor = c
    //intensity = _in
  }



  V3 getShadowRay(V3 Pt){

    V3 shadowRay(position - Pt);
    Dlh = shadowRay.Length();
    shadowRay.Normalize();
    return shadowRay;

  };

  float computShadow(V3 Pt, int index, V3 Ray, int rand_x, int rand_y);
  Color addShadowColor(Color shadColor, float sh){};



private:

  float Dlh;

};






//////////////////////////////////////////////////////////////////////////
////////////////////////////Shapes Classes////////////////////////////////
//////////////////////////////////////////////////////////////////////////

class Shape 
{
public:

  Shader surface;

  Shape(Shader &material){

    surface = material;

  }

  virtual  int intersect(V3 P0, V3 ray){};
  virtual  float getT(){};
  virtual  float getTmax(){};
  virtual  V3 getIntrPt(V3 Pe, float _t, V3 ray){};
  virtual  V3 getNormal(V3 Ph){};
  virtual  Color getShdColor(Light* LightList[], V3 intersectPt, V3 ray){};
  virtual  float getRadius(){};

private:

  V3 current_Normal;

};



/////////////////////////////////SPHERE////////////////////////////////////

class sphere: public Shape
{
public:

  V3 center;
  float radius;

  //===construtor===//
  sphere (V3 pc, float _r, Shader &material):Shape(material)
  {
    center = pc;
    radius = _r;

  }


  //===other functions===//

  virtual  float getRadius(){
   
    return radius;

  }

  //////////////////////////////////////////////////////////////
  //see if intersection happens: return 1 if intersect, 0 if not
  //////////////////////////////////////////////////////////////

  virtual  int intersect(V3 P0,V3 ray) {      
      b = (center-P0)%ray;
      c = (center-P0)%(center-P0)-radius*radius;
      delta = b*b-c;
      
      if(delta>=0){return 1;}
        else {return 0;}
  } 


  ////////////////////////////////////
  //get t num if intersection happens
  ////////////////////////////////////

  virtual  float getT() {
      float t=b-sqrt(delta);
      return t;

  } 


  virtual  float getTmax() {
      float t=b+sqrt(delta);
      return t;

  }


  /////////////////////////////
  //get the intersection point 
  /////////////////////////////  

  virtual  V3 getIntrPt(V3 Pe, float _t, V3 ray) {
      V3 Ph(Pe + _t*ray);
      return Ph;
  } 


  ///////////////////////
  //get the normal at Ph 
  /////////////////////// 
  virtual  V3 getNormal(V3 Ph){
    if(!surface.bmp){
      V3 nh(Ph-center);
      nh = nh/radius;
      current_Normal = nh;
    }
      return current_Normal;
}


  /////////////////////////////////////////////////
  //get the shading color at the intersection point
  /////////////////////////////////////////////////

  virtual  Color getShdColor(Light* LightList[], V3 intersectPt, V3 ray) { 

      Color color_diff = surface.ambColor;
      Color diffColor = surface.diffColor;
  
      V3 nh = (intersectPt-center)/radius; //surface normal

      float u,v;
      float x,y,z;

      if (surface.tx||surface.bmp){
            //===================================//
            //============= get U,V =============//
            //===================================//

            /*------Define Local Coordinate------*/ 
            V3 rx(-1,0,0);  // rotation coordinates
            V3 ry(0,0,1);
            V3 rz(0,-1,0);
            float sx,sy,sz;
            sx = sy = sz = 4*radius;

            /*-----Compute 3D Location X,Y,Z-----*/
            x = (rx%(intersectPt - center))/sx ;
            y = (ry%(intersectPt - center))/sy ;
            z = (rz%(intersectPt - center))/sz ;

            /*------------Comput U,V-------------*/
            float psi = acos(z);
            float yy = y/sin(psi);

            if(yy>1.0){yy = 1.0;}else if(yy<-1.0){yy = -1.0;}

            float theta = acos(yy);
            v = psi/pi;
            u = theta/(2*pi);
            if (x<0) u = 1.0 - u;
      }


      if (surface.tx){
            
          unsigned char dR, dG, dB;

        switch (surface.txType){
          case 1: //============= get Color From TextureMap =============//          
            dR = surface.textureMap[3*(int(v*height)*width + int(u*width))];
            dG = surface.textureMap[3*(int(v*height)*width + int(u*width)) + 1];
            dB = surface.textureMap[3*(int(v*height)*width + int(u*width)) + 2];
            break;

          case 2: //============= get Color From Procedural TextureMap =============//           
            dR = surface.textureMap[3*(int(v*height)*width + int(u*width))];
            dG = surface.textureMap[3*(int(v*height)*width + int(u*width)) + 1];
            dB = surface.textureMap[3*(int(v*height)*width + int(u*width)) + 2];
            break;

          case 3: //============= get Color From SolidTexture =============//

            /*------------Comput U,V-------------*/
            v = z - (int)z;
            u = x - (int)x;
            if (u < 0) u = 1.0 + u;
            if (v < 0) v = 1.0 + v;
           
            dR = surface.textureMap[3*(int(v*height)*width + int(u*width))];
            dG = surface.textureMap[3*(int(v*height)*width + int(u*width)) + 1];
            dB = surface.textureMap[3*(int(v*height)*width + int(u*width)) + 2];
            break;

        }

            Color mapColor(dR,dG,dB);
            diffColor =  mapColor ;

      }


      if (surface.bmp){
          //===================================================//
          //============= get Normal From BumpMap =============//
          //===================================================//
          V3 v0(0,0,0);
          V3 v1(0,0,0);
          float d0, d1;

          /*------Define Local Coordinate------*/ 
          V3 rx(-1,0,0);  // rotation coordinates
          V3 ry(0,0,1);
          V3 rz(0,-1,0);
          float sx,sy,sz;
          sx = sy = sz = radius;

          /*----------get v0----------*/
          float x = sin(2*pi*(u+1.0/width2))*sin(pi*v);
          float y = cos(2*pi*(u+1.0/width2))*sin(pi*v);
          float z = cos(pi*v);
          V3 Pt0(center + sx*x*rx + sy*y*ry + sz*z*rz) ;
          v0 = Pt0 - intersectPt;
          v0.Normalize();
         
          /*----------get v1----------*/
          x = sin(2*pi*u)*sin(pi*(v+1.0/height2));
          y = cos(2*pi*u)*sin(pi*(v+1.0/height2));
          z = cos(pi*(v+1.0/height2));
          V3 Pt1(center + sx*x*rx + sy*y*ry + sz*z*rz) ;
          v1 = Pt1 - intersectPt;
          v1.Normalize();
        
          /*--------get d0, d1--------*/
          d0 = 2*float(surface.bumpMap[3*(int(v*height2)*width2 + int(u*width2))])/255.0 - 1.0;
          d1 = 2*float(surface.bumpMap[3*(int(v*height2)*width2 + int(u*width2)) + 1])/255.0 - 1.0;

          //cout<< d0 << " " << d1<< " ";

          /*------compute normal------*/
          nh = nh + d0*v0 + d1*v1;
          current_Normal = nh;
    
      }

      //===================================//
      //=========Get Diffuse Color=========//
      //===================================//
    
      for(int nL=0; nL<Num_Light; nL++){

         V3 nLh;
         float c = 0;
         float cosL1 = 0;
         float angle = 0;

         float proj_x, proj_y, proj_u, proj_v, max;
         V3 nxx(-1,0,0);
         V3 nyy(0,-1,0);
         V3 n_project(0,0,0);

         switch (LightList[nL]->type){

         case 1: //for pointlight

            nLh = LightList[nL]->position-intersectPt; //the light ray hitting the pt 
            nLh.Normalize();
            cosL1 = nh%nLh ;//angle with the lightsource
            c = (cosL1+a)/(1.0-a);
            c = LightList[nL]->intensity*c;

              if(c<0){c = 0;}
              if(c>1){c = 1;}

            break;

         case 2: //for spotlight

            nLh = LightList[nL]->position-intersectPt; //the light ray hitting the pt 
            nLh.Normalize();
            cosL1 = nh%nLh ;//angle with the lightsource
            angle = - LightList[nL]->direction%nLh ; 
            c = cosL1;
            c = LightList[nL]->intensity*c;

              if(c<0){c = 0;}
              if(c>1){c = 1;}

            if(angle < LightList[nL]->cosAngle ){c=0;}else{c=c;}

            n_project = (nLh%LightList[nL]->direction)*LightList[nL]->direction - nLh;

            proj_x = n_project%nxx;
            proj_y = n_project%nyy;

            max = sqrt(1.0-LightList[nL]->cosAngle*LightList[nL]->cosAngle);

            proj_u = (proj_x/max + 1.0)/2.0;
            proj_v = (proj_y/max + 1.0)/2.0;

            //cout<< proj_u <<" "<<proj_v <<" ";

            break;


         case 3: //for directionalLight

            nLh = -1*LightList[nL]->direction; //the light ray hitting the pt 
            nLh.Normalize();
            cosL1 = nh%nLh ;//angle with the lightsource
            c = (cosL1+a)/(1.0-a);
            c = LightList[nL]->intensity*c;

              if(c<0){c = 0;}
              if(c>1){c = 1;}

            break;


         case 4: //for areaLight
            int SM = (int)(LightList[nL]->Lsx)/(LightList[nL]->LM); //get areaLight sample number
            int SN = (int)(LightList[nL]->Lsy)/(LightList[nL]->LN);

            for(int jj=0; jj<SN; jj++)
              for(int ii=0; ii<SM; ii++){
                nLh = (LightList[nL]->position + (float(ii)/float(SM))*(LightList[nL]->Lsx)*(LightList[nL]->n0) + (float(jj)/float(SN))*(LightList[nL]->Lsy)*LightList[nL]->getn1())-intersectPt; //the light ray hitting the pt 
                nLh.Normalize();
                cosL1 = nh%nLh;//angle with the lightsource
                c = c + cosL1;
              }

            c = LightList[nL]->intensity*c;

              if(c<0){c = 0;}
              if(c>1){c = 1;}

            break;

         }


        //mixed with Light Color

        unsigned char r_new = diffColor.r*(pixmap2[3*(int(proj_v*height2)*width2 + int(proj_u*width2))]/255.0);
        unsigned char g_new = diffColor.r*(pixmap2[3*(int(proj_v*height2)*width2 + int(proj_u*width2))+1]/255.0);
        unsigned char b_new = diffColor.r*(pixmap2[3*(int(proj_v*height2)*width2 + int(proj_u*width2))+2]/255.0);



        //unsigned char r_new = diffColor.r*(LightList[nL]->Lcolor.r/255.0);
        //unsigned char g_new = diffColor.g*(LightList[nL]->Lcolor.g/255.0);
        //unsigned char b_new = diffColor.b*(LightList[nL]->Lcolor.b/255.0);
        Color finalDiffColor(r_new, g_new, b_new);

        color_diff = finalDiffColor.mixColor(color_diff, c);
        // color_diff = diffColor.mixColor(color_diff, c);

      }


         Color color_final = color_diff;




        //====================================//
        //=========Get Sepcular Color=========//
        //====================================//

      for(int nL=0; nL<Num_Light; nL++){

         float s = 0; 
         float a = 6.0;
         float s0 = 0.5;
         float det = 0.48;

         V3 nLh;
         V3 r(ray-2*(ray%nh)*nh);
         r.Normalize();

         switch (LightList[nL]->type){

         case 1: //for pointlight

            nLh = LightList[nL]->position-intersectPt; //the light ray hitting the pt 
            nLh.Normalize();

            s = nLh%r; 
            s = (s-s0)/det;
            s = LightList[nL]->intensity*s;

            if(s<0){s=0;}
              else if(s>1){s=1;}
                 else {s=pow(s,a);}

            break;


         case 2: //for spotlight

            nLh = LightList[nL]->position-intersectPt; //the light ray hitting the pt 
            nLh.Normalize();

            s = nLh%r; 
            s=(s-s0)/det;
            s = LightList[nL]->intensity*s;

            if(s<0){s=0;}
              else if(s>1){s=1;}
                 else {s=pow(s,a);}

            break;



         case 3: //for directionallight

            nLh = -1*LightList[nL]->direction; //the light ray hitting the pt 
            nLh.Normalize();

            s = nLh%r; 
            s = (s-s0)/det;
            s = LightList[nL]->intensity*s;

            if(s<0){s=0;}
              else if(s>1){s=1;}
                 else {s=pow(s,a);}

            break;


         case 4: //for arealight

            float s1=0;
            int SM = (int)(LightList[nL]->Lsx)/(LightList[nL]->LM); //get areaLight sample number
            int SN = (int)(LightList[nL]->Lsy)/(LightList[nL]->LN);

            for(int jj=0; jj<SN; jj++)
              for(int ii=0; ii<SM; ii++){
               nLh = (LightList[nL]->position + (float(ii)/float(SM))*(LightList[nL]->Lsx)*(LightList[nL]->n0) + (float(jj)/float(SN))*(LightList[nL]->Lsy)*LightList[nL]->getn1())-intersectPt; //the light ray hitting the pt 
               nLh.Normalize();
               s = s + nLh%r;
            }

            s = LightList[nL]->intensity*s;
            s=(s-s0)/det;
               if(s<0){s=0;}
                 else if(s>1){s=1;}
                    else {s=pow(s,a);}

            break;
         }

           color_final = surface.specColor.mixColor(color_final, s);
       }

           return color_final;
  } 

private:

  float delta,b,c;
  V3 current_Normal;

};






///////////////////////////////////PLANE//////////////////////////////////////

class plane: public Shape
{
public:

  V3 center;
  V3 normal;

  //===construtor===//
  plane (V3 pc, V3 n, Shader &material):Shape(material)
  {
    center = pc;
    normal = n;

  }

  //===other functions===//

  virtual  float getRadius(){}


  ///////////////////////////////////////////////////////////////
  //see if intersection happens: return 1 if intersect, 0 if not
  ///////////////////////////////////////////////////////////////

  virtual  int intersect(V3 P0, V3 ray) {      
      t = -(normal%(P0-center))/(normal%ray);      
      if(t>=0){return 1;}
        else {return 0;}
  } 


  ///////////////////////////////////
  //get t num if intersection happens
  ///////////////////////////////////

  virtual  float getT() {
      return t;

  } 


  ////////////////////////////
  //get the intersection point
  ////////////////////////////
   
  virtual  V3 getIntrPt(V3 Pe, float _t, V3 ray) {
      V3 Ph(Pe + _t*ray);
      return Ph;
  } 


  virtual  V3 getNormal(V3 Ph){
      return normal;
  }

  /////////////////////////////////////////////////
  //get the shading color at the intersection point
  /////////////////////////////////////////////////

  virtual  Color getShdColor(Light* LightList[], V3 intersectPt, V3 ray) {

   Color color_diff = surface.ambColor;
   Color diffColor = surface.diffColor;
   float u,v;

    if (surface.tx||surface.bmp){

          //===================================//
          //============= get U,V =============//
          //===================================//

          /*------Define Local Coordinate------*/ 
          V3 rx(-1,0,0);
          V3 ry(0,1,0);
          V3 rz(normal);
          float sx,sy,sz;
          sx = sy = sz = 60.0;

          /*-----Compute 3D Location X,Y,Z-----*/
          float x = (rx%(intersectPt - center))/sx ;
          float y = (ry%(intersectPt - center))/sy ;

          /*------------Comput U,V-------------*/
          v = y - (int)y;
          u = x - (int)x;
          if (u < 0) u = 1.0 + u;
          if (v < 0) v = 1.0 + v;
    
    }



    if (surface.tx){

          unsigned char dR, dG, dB;

        switch (surface.txType){
          case 1: //============= get Color From TextureMap =============//         
            dR = surface.textureMap[3*(int(v*height)*width + int(u*width))];
            dG = surface.textureMap[3*(int(v*height)*width + int(u*width)) + 1];
            dB = surface.textureMap[3*(int(v*height)*width + int(u*width)) + 2];
            break;

          case 2: //============= get Color From Procedural TextureMap =============//
            //int ni = 5;
            int du, dv;
            if(sin(2*pi*2*u)>=0){du = 1;}else{du = -1;}
            if(sin(2*pi*2*v)>=0){dv = 1;}else{dv = -1;}
        
            if(du == dv){          
                dR = 160; dG = 160; dB = 50;
            }else{
                dR = 250; dG = 250; dB = 250;
            }
            break;

          case 3: //============= get Color From SolidTexture =============//          
            dR = surface.textureMap[3*(int(v*height)*width + int(u*width))];
            dG = surface.textureMap[3*(int(v*height)*width + int(u*width)) + 1];
            dB = surface.textureMap[3*(int(v*height)*width + int(u*width)) + 2];
            break;

        }

            Color mapColor(dR,dG,dB);
            diffColor =  mapColor ;

    }



      //===================================//
      //=========Get Diffuse Color=========//
      //===================================//

   for(int nL=0; nL<Num_Light; nL++){

      V3 nLh;
      float c;
      float cosL1;
      float angle;


         float proj_x, proj_y, proj_u, proj_v, max;
         V3 nxx(-1,0,0);
         V3 nyy(0,-1,0);
         V3 n_project(0,0,0);

      switch (LightList[nL]->type){

      case 1: //for pointlight

         nLh = LightList[nL]->position-intersectPt; //the light ray hitting the pt 
         nLh.Normalize();
         cosL1 = normal%nLh ;//angle with the lightsource
         if(cosL1<0){cosL1 = 0;}
         c = cosL1;
         c = LightList[nL]->intensity*c;

         break;

      case 2: //for spotlight

         nLh = LightList[nL]->position-intersectPt; //the light ray hitting the pt 
         nLh.Normalize();
         cosL1 = normal%nLh ;//angle with the lightsource
         if(cosL1<0){cosL1 = 0;}
         c = cosL1;

         angle = - LightList[nL]->direction%nLh ; 
         if(angle < LightList[nL]->cosAngle ){c=0;}else{c=c;}

         c = LightList[nL]->intensity*c;


            n_project = (nLh%LightList[nL]->direction)*LightList[nL]->direction - nLh;
            

            proj_x = n_project%nxx;
            proj_y = n_project%nyy;

            max = sqrt(1.0-LightList[nL]->cosAngle*LightList[nL]->cosAngle);

            proj_u = (proj_x/max + 1.0)/2.0;
            proj_v = (proj_y/max + 1.0)/2.0;

         break;


      case 3: //for directionallight

         nLh = -1*LightList[nL]->direction; //the light ray hitting the pt 
         nLh.Normalize();
         cosL1 = normal%nLh ;//angle with the lightsource
         if(cosL1<0){cosL1 = 0;}
         c = cosL1;
         c = LightList[nL]->intensity*c;

         break;


      case 4: //for arealight

         int SM = (int)(LightList[nL]->Lsx)/(LightList[nL]->LM); //get areaLight sample number
         int SN = (int)(LightList[nL]->Lsy)/(LightList[nL]->LN);

         for(int jj=0; jj<SN; jj++)
           for(int ii=0; ii<SM; ii++){
             nLh = (LightList[nL]->position + (float(ii)/float(SM))*(LightList[nL]->Lsx)*(LightList[nL]->n0) + (float(jj)/float(SN))*(LightList[nL]->Lsy)*LightList[nL]->getn1())-intersectPt; //the light ray hitting the pt 
             nLh.Normalize();
             cosL1 = normal%nLh;//angle with the lightsource
             c = c + cosL1;
           }

         c = LightList[nL]->intensity*c;

           if(c<0){c = 0;}
           if(c>1){c = 1;}

         break;

      }

        //mixed with Light Color

        unsigned char r_new = diffColor.r*(pixmap2[3*(int(proj_v*height2)*width2 + int(proj_u*width2))]/255.0);
        unsigned char g_new = diffColor.r*(pixmap2[3*(int(proj_v*height2)*width2 + int(proj_u*width2))+1]/255.0);
        unsigned char b_new = diffColor.r*(pixmap2[3*(int(proj_v*height2)*width2 + int(proj_u*width2))+2]/255.0);

        //unsigned char r_new = diffColor.r*(LightList[nL]->Lcolor.r/255.0);
        //unsigned char g_new = diffColor.g*(LightList[nL]->Lcolor.g/255.0);
        //unsigned char b_new = diffColor.b*(LightList[nL]->Lcolor.b/255.0);
        Color finalDiffColor(r_new, g_new, b_new);

        color_diff = finalDiffColor.mixColor(color_diff, c);
   }


      return color_diff;
} 


private:

  float t;
  V3 current_Normal;
  
};



Shape **GeoList;




///////////////////////////////////////////////////////////////////////////////


float pointLight::computShadow(V3 Pt, int index, V3 Ray, int rand_x, int rand_y){

           float sh;
           float L=0;
           float tmin2,tmax2;
           int time2 = 0;
           float soft_a = 10.0; 

           for(int ng=0; ng<Num; ng++){
            if(ng != index){  
             if (GeoList[ng]->intersect(Pt,Ray)==1){
               time2++;
               tmin2 = GeoList[ng]->getT();
               tmax2 = GeoList[ng]->getTmax();

               float r = GeoList[ng]->getRadius();

               if(tmin2>t_shadowBias && tmin2<Dlh){
               float d = (1.6*(tmax2-tmin2)-0.6)/r;
               L = L + pow(d,soft_a);}
              }} 
           }


           if (time2!=0){
           sh = L;
           if (sh>1){sh = 1.0;}
           } else {sh = 0;}

           return sh;

  };




float directionalLight::computShadow(V3 Pt, int index, V3 Ray, int rand_x, int rand_y){

           float sh=0;
           float L=0;
           int time2 = 0;
           Ray.Normalize();

           for(int ng=0; ng<Num; ng++){
            if(ng != index){  
             if (GeoList[ng]->intersect(Pt,Ray)==1){
               float tmin2 = GeoList[ng]->getT();
               if(tmin2>t_shadowBias){time2++;}
              }
             } 
           }

           if (time2!=0){sh = 1.0;}
           return sh;

  };







float spotLight::computShadow(V3 Pt, int index, V3 Ray,int rand_x, int rand_y){

           float sh;
           float L=0;
           float tmin2,tmax2;
           int time2 = 0;

           for(int ng=0; ng<Num; ng++){
            if(ng != index){  
             if (GeoList[ng]->intersect(Pt,Ray)==1){
               time2++;
               tmin2 = GeoList[ng]->getT();
               tmax2 = GeoList[ng]->getTmax();

               if(tmin2>t_shadowBias && tmin2<Dlh){
               L = L+(tmax2-tmin2);}
              }} 
           }


           if (time2!=0){
           sh = 4*(L/Dlh);
           if (sh>1){sh = 1.0;}
           } else {sh = 0;}

           return sh;

  };



float areaLight::computShadow(V3 Pt, int index, V3 Ray, int rand_x, int rand_y){

           float sh;
           float L=0;
           float tmin2,tmax2;

           V3 ny(direction^n0); 
           ny.Normalize();

           int SM = (int)Lsx/LM; //get areaLight sample number
           int SN = (int)Lsy/LN;

           int time2 = 0;

           float radm = 1.0*rand()/RAND_MAX;
           float radn = 1.0*rand()/RAND_MAX;

           V3 LightPt(position + ((rand_x%M)*LM + LM*radm)*n0 + ((rand_y%N)*LN + radn*LN)*ny);

           V3 shadowRay(LightPt - Pt);
           Dlh = shadowRay.Length();
           shadowRay.Normalize();

           for(int ng=0; ng<Num; ng++){
            if(ng != index){  
             if (GeoList[ng]->intersect(Pt,shadowRay)==1){
               tmin2 = GeoList[ng]->getT();
               if(tmin2>t_shadowBias && tmin2<Dlh){time2++;}
              }
             } 
           }

           if (time2!=0){sh = 1.0;}



           return sh;

  };



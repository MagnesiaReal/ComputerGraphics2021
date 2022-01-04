#ifndef LIGHTING
#define LIGHTING
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include "animation_structs.h"

// atenuation factor
double att (int d, Light l){
  return (double)(1/(l.a*d*d + l.b*d + l.c));
}

double attp(int d) {
  return (double)(1/(0.000189*d*d + 0.000018*d + -0.08));
}

double dist(struct xyz v1, struct xyz v2) {
  double xsq = (v1.x - v2.x)*(v1.x - v2.x);
  double ysq = (v1.y - v2.y)*(v1.y - v2.y);
  double zsq = (v1.z - v2.z)*(v1.z - v2.z);
  return sqrt((double)sqrt(xsq + ysq + zsq));
}



bool is_light(struct xyz normal, struct xyz vect){
 
  double nh = sqrt((double)(normal.x*normal.x + normal.z*normal.z));
  double vh = sqrt ((double)(vect.x*vect.x + vect.z*vect.z));

  double nhy = sqrt((double)(normal.y*normal.y + normal.z*normal.z));
  double vhy = sqrt ((double)(vect.y*vect.y + vect.z*vect.z));

  double nxz = (double)(normal.x/nh);
  double vxz = (double)(vect.x/vh);
  double nyz = (double)(normal.y/nhy);
  double vyz = (double)(vect.y/vhy);
  printf("nxz = %lf  vxz = %lf\n", nxz, vxz);
  
  if(abs(nxz + vxz) < 1 || abs(nyz + vyz) < 1) return false;
  else return true;
}

void ilumination(struct face* faces, unsigned int size, Light l, double ambient) {
  
  struct xyz lp = {l.x, l.y, l.z};

  for (unsigned int f = 0; f < size; ++f) {
  
    if(faces[f].draw) {
      
      if (faces[f].la != NULL) {
        free(faces[f].la);
        faces[f].la = NULL;
      }
      faces[f].la = (Matrix_light*)calloc(faces[f].fullfill_size, sizeof(Matrix_light));
      
      
      struct xyz* array = faces[f].fullfill_px;
      
      for (unsigned int i = 0; i < faces[f].fullfill_size; ++i) {
        
        faces[f].la[i].R = faces[f].R*ambient;
        faces[f].la[i].G = faces[f].G*ambient;
        faces[f].la[i].B = faces[f].B*ambient;

        // break if
        //struct xyz vect_dir = {
           //array[i].x -lp.x,
           //array[i].y -lp.y,
           //array[i].z -lp.z
        //};
        //if (is_light(faces[f].normal_vector, vect_dir)) break;


        double at = att(dist(array[i], lp), l);
        unsigned char lr = l.R*at; 
        unsigned char lg = l.G*at; 
        unsigned char lb = l.B*at; 
        // attenuation factor
        if(lr <= faces[f].R) {
          
          if (lr > (int)(faces[f].R*ambient)) faces[f].la[i].R = lr;

          //printf("att=%lf lr=%d ambient=%i \n",at, lr, (int)(faces[f].R*ambient));
        }
        else faces[f].la[i].R = faces[f].R;

        if(lg <= faces[f].G) {
          if (lg > (int)(faces[f].G*ambient)) faces[f].la[i].G = lg;
          
          //printf("att=%lf lg=%d ambient=%i \n",at, lg, (int)(faces[f].G*ambient));
        }
        else faces[f].la[i].G = faces[f].G;

        if(lb <= faces[f].B) {
          if (lb > (int)(faces[f].B*ambient)) faces[f].la[i].B = lb;
          
          //printf("att=%lf lb=%d ambient=%i \n",at, lb, (int)(faces[f].B*ambient));
        }
        else faces[f].la[i].B = faces[f].B;
      }
    }    
  }
}


void iluminationp(struct face* faces, unsigned int size, int f, double ambient) {
  
  for (unsigned int f = 0; f < size; ++f) {
    if(faces[f].draw){
      for (unsigned int i = 0; i < faces[f].fullfill_size; ++i){
        double at = attp(faces[f].fullfill_px[i].z - f); 
        unsigned char lr = faces[f].la[i].R*at; 
        unsigned char lg = faces[f].la[i].G*at; 
        unsigned char lb = faces[f].la[i].B*at;

        if (at >= 1.0000000) {
          printf("at overflow = %lf\n", at);
        }
        if(lr <= faces[f].R) {
          
          if (lr > (int)(faces[f].R*ambient)) faces[f].la[i].R = lr;

          //printf("att=%lf lr=%d ambient=%i \n",at, lr, (int)(faces[f].R*ambient));
        }
        else faces[f].la[i].R = faces[f].R;

        if(lg <= faces[f].G) {
          if (lg > (int)(faces[f].G*ambient)) faces[f].la[i].G = lg;
          
          //printf("att=%lf lg=%d ambient=%i \n",at, lg, (int)(faces[f].G*ambient));
        }
        else faces[f].la[i].G = faces[f].G;

        if(lb <= faces[f].B) {
          if (lb > (int)(faces[f].B*ambient)) faces[f].la[i].B = lb;
          
          //printf("att=%lf lb=%d ambient=%i \n",at, lb, (int)(faces[f].B*ambient));
        }
        else faces[f].la[i].B = faces[f].B;

      }
      
      

    }
  }
}

#endif

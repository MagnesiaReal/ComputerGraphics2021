#ifndef FACE_HIDDING
#define FACE_HIDDING
#include "animation_structs.h"
#include <math.h>
#include <stdio.h>


struct xyz cross_product(struct xyz vec1, struct xyz vec2) {
  struct xyz abac = {
    vec1.y*vec2.z - vec1.z*vec2.y,
    vec1.z*vec2.x - vec1.x*vec2.z,
    vec1.x*vec2.y - vec1.y*vec2.x,
  };
  return abac;
}


// you are working with x pointing to opositive direction
// so, cross product give you real direction of Z axis
// in conclusion: you get the positive Z like negative Z.
void save_normal_vector(struct face *faces, struct edge *edges, double **vertex, unsigned int faces_size) {
  
  for (unsigned int f = 0; f < faces_size; ++f) {
    struct xyz ab = {
      (int)(vertex[edges[faces[f].e1-1].v2-1][0] - vertex[edges[faces[f].e1-1].v1-1][0]),
      (int)(vertex[edges[faces[f].e1-1].v2-1][1] - vertex[edges[faces[f].e1-1].v1-1][1]),
      (int)(vertex[edges[faces[f].e1-1].v2-1][2] - vertex[edges[faces[f].e1-1].v1-1][2])  
    };
    struct xyz ac = {
     (int)(vertex[edges[faces[f].e2-1].v2-1][0] - vertex[edges[faces[f].e1-1].v1-1][0]),
     (int)(vertex[edges[faces[f].e2-1].v2-1][1] - vertex[edges[faces[f].e1-1].v1-1][1]),
     (int)(vertex[edges[faces[f].e2-1].v2-1][2] - vertex[edges[faces[f].e1-1].v1-1][2])
    };
    //printf("ab : {x : %i, y : %i, z : %i}\n", ab.x, ab.y , ab.z);
    //printf("ac : {x : %i, y : %i, z : %i}\n", ac.x, ac.y , ac.z);
    struct xyz abac = cross_product(ab, ac);
    //printf("faces[%i] = {abac : {x : %i, y : %i, z : %i}}\n",f, abac.x, abac.y , abac.z);
    faces[f].normal_vector = abac;
    // calculating if draw or not CORRECT FORM
    //double xz_angle = (abac.z != 0) ? cos((2*M_PI)*abac.x/abac.z) : 0, 
           //yz_angle = (abac.z != 0) ? cos((2*M_PI)*abac.y/abac.z) : 0;
    
    //faces[f].draw = (xz_angle >= 0 && yz_angle >= 0);
    //printf("face %i: {abac: {x: %i, y: %i, z: %i}, acos xz: %lf, acos yz: %lf, draw: %s}\n", f, abac.x, abac.y, abac.z, xz_angle, yz_angle, (faces[f].draw)?"true":"false");
    
    // calculating if draw or not DIRTY FORM
    faces[f].draw = !(abac.z <= 0); //this is negate cause it's pointing to real direction Z.
  }
}
#endif

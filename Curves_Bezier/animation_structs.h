#ifndef ANIMATION_STRUCTS
#define ANIMATION_STRUCTS
#include <stdbool.h>

typedef struct light {
  // position
  int x, y ,z;
  // parameters
  double a, b, c;
  // color
  unsigned char R, G, B;
} Light;

typedef struct matrix_light {
  unsigned char R, G, B;
} Matrix_light;



struct xyz {
  int x, y, z;
};

struct edge {
  unsigned int v1, v2;
};

struct face {
  unsigned int e1;
  unsigned int e2;
  unsigned int e3;

  unsigned int px_size, fullfill_size;
  struct xyz *px_array, *fullfill_px, normal_vector;
  bool draw;

  // color const 
  unsigned char R, G, B;
  // color affected by lightsource
  
  Matrix_light *la;
  // diffuse factor
  double dr, dg, db; // no greatter than one
};

#endif

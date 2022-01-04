#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <float.h>
#include <unistdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include "face_hidding.h"

struct xyzd{
  double x,y,z;
};

// FOR PPM /////////
Matrix_light **matrix = NULL;
////////////////////

////////////////////
double max_x = -DBL_MAX, min_x = DBL_MAX, max_y = -DBL_MAX, min_y = DBL_MAX, max_z = -DBL_MAX,
       min_z = DBL_MAX, minpos_x = DBL_MAX, minpos_y = DBL_MAX, minpos_z = DBL_MAX;

// ##############################################################
void Bresenham_draw_line (int x1, int y1, int x2, int y2)
{
  int dy = y2 - y1, dx = x2 - x1, diagy, diagx;

  if (dy >= 0)
    diagy = 1;			// north diagonal
  else
  {				// south diagonal
    dy = -dy;
    diagy = -1;
  }

  if (dx >= 0)
    diagx = 1;			// east diagonal
  else
  {				// west diagonal
    dx = -dx;
    diagx = -1;
  }
  int x, y;
  int *a = &x, *b = &y;
  if (dy > dx) {
    int aux = dx;
    dx = dy;
    dy = aux;

    aux = x1;
    x1 = y1;
    y1 = aux;

    aux = x2;
    x2 = y2;
    y2 = aux;

    aux = diagx;
    diagx = diagy;
    diagy = aux;
    a = &y;
    b = &x;
  }
  x = x1;
  y = y1;
  int d = 2 * dy - dx, dup1 = 2 * dy, dup2 = 2 * (dy - dx);
  
  int red=255, green=255 , blue=255;

  while (x != x2)
  {

    matrix[*b][*a].R = red;
    matrix[*b][*a].G = green;
    matrix[*b][*a].B = blue;
    red = 12;
    green =90;
    blue = 199;

    if (d >= 0)
    {
      x += diagx;
      y += diagy;
      d += dup2;
    }
    else
    {
      x += diagx;
      d += dup1;
    }
  } 
}

void Bresenham (struct xyz *ar, int size)
{
  for (int i = 0; i + 1 < size; ++i)
  {
    Bresenham_draw_line(ar[i].x, ar[i].y, ar[i + 1].x, ar[i + 1].y); 
  }
}


void raster() {
  printf("Rasterising\n");
  // raster for edges"q

}

void matrix_multiplication (float A[4][4], float B[4][4], float C[4][4])
{
  // init 0
  for (int i = 0; i < 4; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      C[i][j] = 0;
    }
  }
  // multiplication
  for (int i = 0; i < 4; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      for (int k = 0; k < 4; ++k)
      { 
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}


void Save_ppm (char *filename, int size_x, int size_y)
{
  FILE* ppm_file = fopen (filename, "w");
  if (ppm_file == NULL)
  {
    perror ("No se pudo crear el archivo : ");
    exit (-1);
  }

  fprintf (ppm_file, "P3\n%i %i\n255", (int) (size_x), (int) (size_y));
  for (int i = size_y - 1; i >= 0; --i)
  {
    fputs ("\n", ppm_file);
    for (int j = 0; j < size_x; ++j)
    {
      fprintf (ppm_file, "%d %d %d ", matrix[i][j].R, matrix[i][j].G, matrix[i][j].B);
    }
  }
  fclose (ppm_file);
}


void matrix_C(double A[4], double B[4][3], double C[3]){
  for (int i = 0; i < 3; ++i) {
    for (int k = 0; k < 4; ++k) {
      C[i] += A[k] * B[k][i];
    }
  }
}

void matrix_U(double A[4], struct xyzd B[4][4], struct xyzd C[4]) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      C[i].x += A[j] * B[i][j].x;
      C[i].y += A[j] * B[i][j].y;
      C[i].z += A[j] * B[i][j].z;
    }
  }
}

void matrix_V(double A[4], struct xyzd B[4], struct xyzd* C) {
  for (int i = 0; i < 4; ++i) {
    C->x += A[i] * B[i].x;
    C->y += A[i] * B[i].y;
    C->z += A[i] * B[i].z;

  } 
}

void matrix_multiplication_4x4_4x3 (double A[4][4], double B[4][3], double C[4][3]) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 4; ++k) {
        C[i][j] += A[i][k] * B[k][j];
      }  
    }
  }
}


void  new_vertexes(double P[4][4], struct xyzd **array, int size) {
  
  for (int vh = 0; vh < size; ++vh) 
    for (int v = 0; v < size; ++v) {
      double new_vertex[4] = {0};
      for (int i = 0; i < 4; ++i) {
        new_vertex[i] += P[i][0] * array[vh][v].x;
        new_vertex[i] += P[i][1] * array[vh][v].y;
        new_vertex[i] += P[i][2] * array[vh][v].z;
        new_vertex[i] += P[i][3] * 1;
      }
      
      array[vh][v].x = new_vertex[0] / new_vertex[3];
      array[vh][v].y = new_vertex[1] / new_vertex[3];
      array[vh][v].z = new_vertex[2]; // no perform operation in Z. if occurs z always be f.

    }
}

int main (int argc, char *argv[])
{
  
  if (argc < 48)
  {
    printf("\033[1m\033[39mSurface needs 16 points(4 points and 12 bezier points) with X Y Z values, total 48 parameters:\n\n");
    printf("./surface \\\n\
\033[1m\033[31m<P1x> <P1y> <P1z> \033[0m\033[32m<B1x> <B1y> <B1z> <B2x> <B2y> <B2z> \033[1m\033[31m<P2x> <P2y> <P2z> \\\n\
\033[0m\033[32m<B3x> <B3y> <B3z> \033[0m\033[33m<B5x> <B5y> <B5z> <B7x> <B7y> <B7z> \033[0m\033[32m<B9x> <B9y> <B9z> \\\n\
\033[0m\033[32m<B4x> <B4y> <B4z> \033[0m\033[33m<B6x> <B6y> <B6z> <B8x> <B8y> <B8z> \033[0m\033[32m<B10x> <B10y> <B10z> \\\n\
\033[1m\033[31m<P3x> <P3y> <P3z> \033[0m\033[32m<B11x> <B11y> <B11z> <B12x> <B12y> <B12z> \033[1m\033[31m<P4x> <P4y> <P4z>\033[36m\n\n");
    puts("this program uses 40 fixed lines for each curve.");
    exit (EXIT_FAILURE);
  }

  // Geometrical matrix
  char *string;
  struct xyzd G[4][4] = {
    {
      {strtod(argv[ 1],&string),strtod( argv[ 2],&string),strtod( argv[ 3],&string)}, 
      {strtod(argv[ 4],&string),strtod( argv[ 5],&string),strtod( argv[ 6],&string)}, 
      {strtod(argv[ 7],&string),strtod( argv[ 8],&string),strtod( argv[ 9],&string)}, 
      {strtod(argv[10],&string),strtod( argv[11],&string),strtod( argv[12],&string)}  
    },
    {
      {strtod(argv[13],&string),strtod( argv[14],&string),strtod( argv[15],&string)}, 
      {strtod(argv[16],&string),strtod( argv[17],&string),strtod( argv[18],&string)}, 
      {strtod(argv[19],&string),strtod( argv[20],&string),strtod( argv[21],&string)}, 
      {strtod(argv[22],&string),strtod( argv[23],&string),strtod( argv[24],&string)}
    },
    {
      {strtod(argv[25],&string),strtod( argv[26],&string),strtod( argv[27],&string)}, 
      {strtod(argv[28],&string),strtod( argv[29],&string),strtod( argv[30],&string)}, 
      {strtod(argv[31],&string),strtod( argv[32],&string),strtod( argv[33],&string)}, 
      {strtod(argv[34],&string),strtod( argv[35],&string),strtod( argv[36],&string)}
    },
    {
      {strtod(argv[37],&string),strtod( argv[38],&string),strtod( argv[39],&string)}, 
      {strtod(argv[40],&string),strtod( argv[41],&string),strtod( argv[42],&string)}, 
      {strtod(argv[43],&string),strtod( argv[44],&string),strtod( argv[45],&string)}, 
      {strtod(argv[46],&string),strtod( argv[47],&string),strtod( argv[48],&string)}
    }
  };

  int n_lines = 40;
  int vertexes = n_lines + 1;
 
  struct xyzd **array = malloc(vertexes * sizeof(struct xyzd *));
  for (int i = 0; i < vertexes; ++i) {
    array[i] = calloc(vertexes, sizeof(struct xyzd));
  }
  
  struct xyz **arrayint = malloc(vertexes * sizeof(struct xyz*));
  for (int i = 0; i < vertexes; ++i) {
    arrayint[i] = calloc(vertexes, sizeof(struct xyz));
  }

  double u = 0, du = (double)(1.f/n_lines);
  for (int vertexh = 0; vertexh < vertexes; ++vertexh) {
    double U[4] = {u*u*u, 3*u*u*(1.f-u), 3*u*(1.f-u)*(1.f-u), (1.f-u)*(1.f-u)*(1.f-u)};
    double v = 0, dv = (double)(1.f/n_lines);
    for (int vertexv = 0; vertexv < vertexes; ++vertexv) {
      double V[4] = {v*v*v, 3*v*v*(1.f-v), 3*v*(1.f-v)*(1.f-v), (1.f-v)*(1.f-v)*(1.f-v)}; 
      struct xyzd C[4] = {0};
      matrix_U(U, G, C);
      struct xyzd S = {0};
      matrix_V(V, C, &S);
      array[vertexh][vertexv].x = S.x;
      array[vertexh][vertexv].y = S.y;
      array[vertexh][vertexv].z = S.z;

      if (max_x < S.x) max_x = S.x;
      if (min_x > S.x) min_x = S.x;
      if (S.x >= 0 && minpos_x > S.x) minpos_x = S.x;

      if (max_y < S.y) max_y = S.y;
      if (min_y > S.y) min_y = S.y; 
      if (S.y >= 0 && minpos_y > S.y) minpos_y = S.y;

      if (max_z < S.z) max_z = S.z;
      if (min_z > S.z) min_z = S.z; 
      if (S.z >= 0 && minpos_z > S.z) minpos_z = S.z;

      // increment
      printf("u=%lf v=%lf  x=%lf y=%lf z=%lf\n",u,v,S.x,S.y,S.z);
      v += dv;
      
    }
    // increment
    u += du;
  }
  double f = 0.8f*minpos_z;
  // perpective projection  matrix
  double Pp[4][4] = {
    {f, 0, 0, 0},
    {0, f, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 1, 0},
  };

  new_vertexes(Pp, array, vertexes);

  max_x = -DBL_MAX, min_x = DBL_MAX, max_y = -DBL_MAX, min_y = DBL_MAX, max_z = -DBL_MAX,
  min_z = DBL_MAX, minpos_x = DBL_MAX, minpos_y = DBL_MAX, minpos_z = DBL_MAX;
  for (int vertexh = 0; vertexh < vertexes; ++vertexh) 
    for (int vertexv = 0; vertexv < vertexes; ++vertexv) {
      double C[3] = {array[vertexh][vertexv].x, array[vertexh][vertexv].y,array[vertexh][vertexv].z};

      if (max_x < C[0]) max_x = C[0];
      if (min_x > C[0]) min_x = C[0];
      if (C[0] >= 0 && minpos_x > C[0]) minpos_x = C[0];

      if (max_y < C[1]) max_y = C[1];
      if (min_y > C[1]) min_y = C[1]; 
      if (C[1] >= 0 && minpos_y > C[1]) minpos_y = C[1];

      if (max_z < C[2]) max_z = C[2];
      if (min_z > C[2]) min_z = C[2]; 
      if (C[2] >= 0 && minpos_z > C[2]) minpos_z = C[2];

      //printf("x=%lf y=%lf z=%lf\n", C[0], C[1], C[2]);
    }

  int mul;
  if (abs((int)(max_x - min_x)) < 25 || abs((int)(max_y - min_y)) < 25) mul = 100; 
  else if (abs((int)(max_x - min_x)) < 250 || abs((int)(max_y - min_y)) < 250) mul = 10;
  else mul = 1;
  
  printf("multiplier = %i\n\n", mul);

  int mx = (int) ((min_x - minpos_x) * mul);
  int Mx = (int) (max_x * mul);

  int my = (int) ((min_y - minpos_y) * mul);
  int My = (int) (max_y * mul);
  
  for (int vertexh = 0; vertexh < vertexes; ++vertexh) 
    for (int vertexv = 0; vertexv < vertexes; ++vertexv) {
      arrayint[vertexh][vertexv].x = (int)((array[vertexh][vertexv].x - min_x + minpos_x) * mul);
      arrayint[vertexh][vertexv].y = (int)((array[vertexh][vertexv].y - min_y + minpos_y) * mul);
      arrayint[vertexh][vertexv].z = (int)(array[vertexh][vertexv].z * mul);
      printf("x=%i y=%i z=%i\n", arrayint[vertexh][vertexv].x, arrayint[vertexh][vertexv].y, arrayint[vertexh][vertexv].z);
    } 

  printf("matrix [%i][%i]\n", My - my + 2, Mx - mx + 2);
  matrix = malloc((My - my + 2)*sizeof(Matrix_light*));
  for (int i = 0; i < My - my + 2; ++i) {
    matrix[i] = calloc(Mx - mx + 2, sizeof(Matrix_light));   
  }
  
  struct xyz **arrayinth = malloc(vertexes * sizeof(struct xyz*));
  for (int vertex = 0; vertex < vertexes; ++vertex) {
    arrayinth[vertex] = calloc(vertexes, sizeof(struct xyz));
    for (int v = 0; v < vertexes; ++v) {
      arrayinth[vertex][v].x = arrayint[v][vertex].x;
      arrayinth[vertex][v].y = arrayint[v][vertex].y;
      arrayinth[vertex][v].z = arrayint[v][vertex].z;
    }
  }

  for (int vertex = 0; vertex < vertexes; ++vertex) {
    Bresenham(arrayint[vertex], vertexes);
    Bresenham(arrayinth[vertex], vertexes);
  }
  
  // raster();
  Save_ppm("surface.ppm", Mx - mx + 2, My - my + 2);

  return 0;
}

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

void matrix_multiplication_4x4_4x3 (double A[4][4], double B[4][3], double C[4][3]) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 4; ++k) {
        C[i][j] += A[i][k] * B[k][j];
      }  
    }
  }
}


void  new_vertexes(double P[4][4], struct xyzd array[], int size) {
  
  for (int v = 0; v < size; ++v) {
    double new_vertex[4] = {0};
    for (int i = 0; i < 4; ++i) {
      new_vertex[i] += P[i][0] * array[v].x;
      new_vertex[i] += P[i][1] * array[v].y;
      new_vertex[i] += P[i][2] * array[v].z;
      new_vertex[i] += P[i][3] * 1;
    }
    
    array[v].x = new_vertex[0] / new_vertex[3];
    array[v].y = new_vertex[1] / new_vertex[3];
    array[v].z = new_vertex[2]; // no perform operation in Z. if occurs z always be f.

  }
}

int main (int argc, char *argv[])
{
  
  if (argc < 12)
  {
    printf("Hermite curve needs four parameters:\n");
    printf("./hermite <p1x> <p1y> <p1z> <p2x> <p2y> <p2z> <h1x> <h1y> <h1z> <h2x> <h2y> <h2z>\n");
    puts("this program uses 200 fixed lines for the curve");
    exit (EXIT_FAILURE);
  }

  /* Curves
   * P(0) = p_k
   * P(1) = p_{k + 1}
   * P'(0) = Dp_k (derivative)
   * P'(1) = Dp_{k + 1} (derivative)
   * vvvvvvvvvvvvvvvvvvvvvvvvv
   * P(u) = au³ + bu² + cu + d
   * P'(u) = 3au² + 2bu + c + d
   * vvvvvvvvvvvvvvvvvvvvvvvvv
   * [ P(0)]   [0 0 0 1] [a]
   * [ P(1)] = [1 1 1 1] [b]
   * [P'(0)]   [0 0 1 0] [c]
   * [P'(1)]   [3 2 1 0] [d]
   * 
   * Hermite Matrix
   *       [0 0 0 1] 
   * M_h = [1 1 1 1]
   *       [0 0 1 0]
   *       [3 2 1 0]
   *
   *         [ 2 -2  1  1]
   * M_h⁻¹ = [-3  3 -2 -1]
   *         [ 0  0  1  0]
   *         [ 1  0  0  0]
   * */
  double Mh[4][4] = { 
    { 2,-2, 1, 1},
    {-3, 3,-2,-1},
    { 0, 0, 1, 0},
    { 1, 0, 0, 0}
  };
  // Geometrical matrix
  char *string;
  double G[4][3] = {
    {strtod(argv[ 1],&string),strtod( argv[ 2],&string),strtod( argv[ 3],&string)}, // Point 1
    {strtod(argv[ 4],&string),strtod( argv[ 5],&string),strtod( argv[ 6],&string)}, // Point 2 
    {strtod(argv[ 7],&string),strtod( argv[ 8],&string),strtod( argv[ 9],&string)}, // Vector 1 Velocity
    {strtod(argv[10],&string),strtod( argv[11],&string),strtod( argv[12],&string)}  // Vector 2 Velocity
  };

  double A[4][3] = {0};
  
  matrix_multiplication_4x4_4x3(Mh, G, A);
  
  // matrix for coefficients 


  int n_lines = 200;
  struct xyzd *array = calloc(n_lines + 1, sizeof(struct xyzd));
  struct xyz *arrayint = calloc(n_lines + 1, sizeof(struct xyz));

  int vertexes = n_lines + 1;

  double t = 0, dt = (double)(1.f/n_lines);
  for (int vertex = 0; vertex < vertexes; ++vertex) {
    // matrix T
    double T[4] = {pow(t,3), pow(t,2), t, 1};
    // matrix Curve
    double C[3] = {0};
    // matrix multiplication
    matrix_C(T, A, C);
    
    printf("t = %lf\t x=%lf y=%lf z=%lf\n", t, C[0], C[1], C[2]);


    
    array[vertex].x = C[0];
    array[vertex].y = C[1];
    array[vertex].z = C[2];

    if (max_x < C[0]) max_x = C[0];
    if (min_x > C[0]) min_x = C[0];
    if (C[0] >= 0 && minpos_x > C[0]) minpos_x = C[0];

    if (max_y < C[1]) max_y = C[1];
    if (min_y > C[1]) min_y = C[1]; 
    if (C[1] >= 0 && minpos_y > C[1]) minpos_y = C[1];

    if (max_z < C[2]) max_z = C[2];
    if (min_z > C[2]) min_z = C[2]; 
    if (C[2] >= 0 && minpos_z > C[2]) minpos_z = C[2];
    
    t += dt;
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

  for (int vertex = 0; vertex < vertexes; ++vertex) {
    double C[3] = {array[vertex].x, array[vertex].y,array[vertex].z};

    if (max_x < C[0]) max_x = C[0];
    if (min_x > C[0]) min_x = C[0];
    if (C[0] >= 0 && minpos_x > C[0]) minpos_x = C[0];

    if (max_y < C[1]) max_y = C[1];
    if (min_y > C[1]) min_y = C[1]; 
    if (C[1] >= 0 && minpos_y > C[1]) minpos_y = C[1];

    if (max_z < C[2]) max_z = C[2];
    if (min_z > C[2]) min_z = C[2]; 
    if (C[2] >= 0 && minpos_z > C[2]) minpos_z = C[2];

    printf("x=%lf y=%lf z=%lf\n", C[0], C[1], C[2]);
  }

  int mul;
  if (abs((int)(max_x - min_x)) < 25 || abs((int)(max_y - min_y)) < 25) mul = 150; 
  else if (abs((int)(max_x - min_x)) < 250 || abs((int)(max_y - min_y)) < 250) mul = 15;
  else mul = 1;
  
  printf("multiplier = %i\n\n", mul);

  int mx = (int) ((min_x - minpos_x) * mul);
  int Mx = (int) (max_x * mul);

  int my = (int) ((min_y - minpos_y) * mul);
  int My = (int) (max_y * mul);
  
  for (int vertex = 0; vertex < vertexes; ++vertex) {
    arrayint[vertex].x = (int)((array[vertex].x - min_x + minpos_x) * mul);
    arrayint[vertex].y = (int)((array[vertex].y - min_y + minpos_y) * mul);
    arrayint[vertex].z = (int)(array[vertex].z * mul);
    printf("x=%i y=%i z=%i\n", arrayint[vertex].x, arrayint[vertex].y, arrayint[vertex].z);
  } 

  printf("matrix [%i][%i]\n", My - my + 2, Mx - mx + 2);
  matrix = malloc((My - my + 2)*sizeof(Matrix_light*));
  for (int i = 0; i < My - my + 2; ++i) {
    matrix[i] = calloc(Mx - mx + 2, sizeof(Matrix_light));   
  }
  

  Bresenham(arrayint, vertexes);
  // raster();
  Save_ppm("hermite_curve.ppm", Mx - mx + 2, My - my + 2);

  return 0;
}

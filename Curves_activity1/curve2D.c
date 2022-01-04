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

// FOR PPM /////////
Matrix_light **matrix = NULL;
////////////////////

////////////////////
double max_x = -DBL_MAX, min_x = DBL_MAX, max_y = -DBL_MAX, min_y = DBL_MAX;

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
  // raster for edges

}

void matrix_multiplication (float A[4][4], float B[4][4], float C[4][4])
{

  for (int i = 0; i < 4; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      C[i][j] = 0;
    }
  }

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

void new_vertex (float final_matrix[4][4])
{
  /*float nv[4] = { 0, 0, 0, 0 };*/
  /*for (int v = 0; v < size_v; ++v)*/
  /*{*/
    /*// vertex operation*/
    /*for (int i = 0; i < 4; ++i)*/
    /*{*/
      /*for (int k = 0; k < 4; ++k)*/
      /*{*/
        /*nv[i] += final_matrix[i][k] * vertexesd[v][k];*/
      /*}*/
    /*}*/
    /*// save new vertex*/
    /*for (int i = 0; i < 4; ++i)*/
    /*{*/
      /*vertexesd[v][i] = nv[i] / nv[3]; // Transform to no homogenious*/
      /*vertexes[v][i] = nv[i] / nv[3];	// copy that*/
      
      /*nv[i] = 0;*/
    /*}*/
  /*}*/
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

void matrix_C(double A[4], double B[4][2], double C[2]){
  for (int i = 0; i < 2; ++i) {
    for (int k = 0; k < 4; ++k) {
      C[i] += A[k] * B[k][i];
    }
  }
}

int main (int argc, char *argv[])
{
  if (argc < 2)
  {
    printf ("Use:\n./curves <ax> <bx> <cx> <dx> <ay> <by> <cy> <dy> <n-lines>");
    exit (EXIT_FAILURE);
  }
  
  // matrix for coefficients 
  char *string;
  double A[4][2] = {
    {strtod(argv[1],&string),strtod( argv[5],&string)},
    {strtod(argv[2],&string),strtod( argv[6],&string)},
    {strtod(argv[3],&string),strtod( argv[7],&string)},
    {strtod(argv[4],&string),strtod( argv[8],&string)}
  };


  int n_lines = strtod(argv[9],&string);
  struct xyzd {
    double x, y, z;
  } *array = calloc(n_lines + 1, sizeof(struct xyzd));
  struct xyz *arrayint = calloc(n_lines + 1, sizeof(struct xyz));

  int vertexes = n_lines + 1;

  double t = 0, dt = (double)(1.f/n_lines);
  for (int vertex = 0; vertex < vertexes; ++vertex) {
    // matrix T
    double T[4] = {pow(t,3), pow(t,2), t, 1};
    // matrix Curve
    double C[2] = {0};
    // matrix multiplication
    matrix_C(T, A, C);
    
    printf("t = %lf\t x=%lf y=%lf\n", t, C[0], C[1]);
    
    array[vertex].x = C[0];
    array[vertex].y = C[1];
    if (max_x < C[0]) max_x = C[0];
    if (max_y < C[1]) max_y = C[1];
    if (min_x > C[0]) min_x = C[0];
    if (min_y > C[1]) min_y = C[1];
    
    t += dt;
  }
  

  int mul;

  if (abs((int)(max_x - min_x)) < 10 || abs((int)(max_y - min_y)) < 10) mul = 100; 
  else if (abs((int)(max_x - min_x)) < 100 || abs((int)(max_y - min_y)) < 100) mul = 10;
  else mul = 1;
  
  printf("multiplier = %i\n\n", mul);

  int mx = (int) (min_x * mul);
  int Mx = (int) (max_x * mul);

  int my = (int) (min_y * mul);
  int My = (int) (max_y * mul);
  
  for (int vertex = 0; vertex < vertexes; ++vertex) {
    arrayint[vertex].x = (int)((array[vertex].x - min_x) * mul);
    arrayint[vertex].y = (int)((array[vertex].y - min_y) * mul);
  } 

  
  matrix = malloc((My - my + 1)*sizeof(Matrix_light*));
  for (int i = 0; i < My - my + 1; ++i) {
    matrix[i] = calloc(Mx - mx + 1, sizeof(Matrix_light));   
  }
  
  struct xyz *exp = calloc(2, sizeof(struct xyz));
  exp[0].x = 1;
  exp[0].y = 1;
  exp[1].x = 888;
  exp[1].y = 777;
  Bresenham(arrayint, vertexes);
  // raster();
  Save_ppm("curve2D.ppm", Mx - mx + 1, My - my + 1);

  return 0;
}

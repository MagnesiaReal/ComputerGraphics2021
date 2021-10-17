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

#define Ver 2001000
#define Tra 4000000

struct args
{
  int ini, end;
};

char bufferus[128];
int nimg = 0;
int maxx = INT_MIN, maxy = INT_MIN;
FILE *ppm_file;
float **vertexes;
int **tray;
int **matrix;
int size_v = 0, size_tray = 0;
float minusx = 0, minusy = 0, minusz = 0, maxix = -FLT_MAX, maxiy =
  -FLT_MAX, maxiz = -FLT_MAX, miniy = FLT_MAX, minix = FLT_MAX, miniz =
  FLT_MAX;
pthread_t *threads;

// ##############################################################
void
Bresenham_draw_line (int x1, int y1, int x2, int y2)
{
  int dy = y2 - y1, dx = x2 - x1, diagy, diagx, recty, rectx;
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

  if (dx >= dy)
    {				// si es mayor la distancia en x que en y entonces
      recty = 0;
      rectx = diagx;
    }
  else
    {
      rectx = 0;
      recty = diagy;
      int temp = dx;
      dx = dy;
      dy = temp;
    }
  int x = x1, y = y1, d = 2 * dy - dx, dup1 = 2 * dy, dup2 = 2 * (dy - dx);

  while (x != x2)
    {
      //printf("matrix[%i][%i]", x, y);
      //fflush(stdout);
      if ((y >= 0 && x >= 0) && (y <= maxiy && x <= maxix))
	matrix[y][x] = 99;

      if (d >= 0)
	{
	  x += diagx;
	  y += diagy;
	  d += dup2;
	}
      else
	{
	  x += rectx;
	  y += recty;
	  d += dup1;
	}
    }
}

void
Bresenham ()
{
  for (int i = 0; i < Tra; ++i)
    {
      if (tray[i][0] == 0)
	break;
      Bresenham_draw_line (vertexes[tray[i][0] - 1][0],
			   vertexes[tray[i][0] - 1][1],
			   vertexes[tray[i][1] - 1][0],
			   vertexes[tray[i][1] - 1][1]);
      Bresenham_draw_line (vertexes[tray[i][1] - 1][0],
			   vertexes[tray[i][1] - 1][1],
			   vertexes[tray[i][2] - 1][0],
			   vertexes[tray[i][2] - 1][1]);
      Bresenham_draw_line (vertexes[tray[i][2] - 1][0],
			   vertexes[tray[i][2] - 1][1],
			   vertexes[tray[i][0] - 1][0],
			   vertexes[tray[i][0] - 1][1]);
      //printf("Figure %i complete\n", i);
    }
}

void *
Bresenham_thread (void *argsp)
{
  struct args *arguments = (struct args *) argsp;
  int i = arguments->ini, end = arguments->end;
  //printf("ini=%i end=%i\n", i, end);
  while (i < end)
    {
      if (tray[i][0] == 0)
	break;
      Bresenham_draw_line (vertexes[tray[i][0] - 1][0],
			   vertexes[tray[i][0] - 1][1],
			   vertexes[tray[i][1] - 1][0],
			   vertexes[tray[i][1] - 1][1]);
      Bresenham_draw_line (vertexes[tray[i][1] - 1][0],
			   vertexes[tray[i][1] - 1][1],
			   vertexes[tray[i][2] - 1][0],
			   vertexes[tray[i][2] - 1][1]);
      Bresenham_draw_line (vertexes[tray[i][2] - 1][0],
			   vertexes[tray[i][2] - 1][1],
			   vertexes[tray[i][0] - 1][0],
			   vertexes[tray[i][0] - 1][1]);
      //printf("Figure %i complete\n", i);
      i++;
    }
  return 0;
}

// ##################################################3


void
matrix_multiplication (float A[4][4], float B[4][4], float C[4][4])
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

void
new_vertex (float final_matrix[4][4])
{
  float nv[4] = { 0, 0, 0, 0 };
  for (int v = 0; v < size_v; ++v)
    {
      // vertex operation
      for (int i = 0; i < 4; ++i)
	{
	  for (int k = 0; k < 4; ++k)
	    {
	      nv[i] += final_matrix[i][k] * vertexes[v][k];
	    }
	}
      // save new vertex
      for (int i = 0; i < 4; ++i)
	{
	  vertexes[v][i] = nv[i] / nv[3];	// tranform to no homogenious
	  nv[i] = 0;
	}
    }
}


void
Save_ppm (const char *filename)
{
  ppm_file = fopen (filename, "w");
  if (ppm_file == NULL)
    {
      perror ("No se pudo crear el archivo : ");
      exit (-1);
    }
  fprintf (ppm_file, "P3\n%i %i\n255", (int) maxix + 1, (int) maxiy + 1);
  for (int i = maxiy; i >= 0; --i)
    {
      fputs ("\n", ppm_file);
      for (int j = 0; j <= maxix; ++j)
	{
	  fprintf (ppm_file, "%i %i %i ", matrix[i][j], matrix[i][j],
		   matrix[i][j]);
	}
    }
  fclose (ppm_file);
}

int
main (int argc, char *argv[])
{
  FILE *obj_file;
  if (argc < 2)
    {
      printf ("Use:\n./objtoppm <objfile>");
      exit (EXIT_FAILURE);
    }
  obj_file = fopen (argv[1], "r");
  if (obj_file == NULL)
    {
      perror ("No se pudo abrir el archivo ");
      exit (-1);
    }
  // Building arrays of arrays for vertexes (x, y) and trays (v1, v2, v3)
  vertexes = (float **) malloc (Ver * sizeof (float *));
  for (int itr = 0; itr < Ver; itr++)
    {
      vertexes[itr] = (float *) malloc (4 * sizeof (float));
      vertexes[itr][3] = 1;
    }

  tray = (int **) malloc (Tra * sizeof (int *));
  for (int itr = 0; itr < Tra; itr++)
    {
      tray[itr] = (int *) malloc (3 * sizeof (int));
    }

  char line[128];
  int i = 0, j = 0;
  float sef = 60.f;
  while (fgets (line, sizeof (line), obj_file) != NULL)
    {
      //fputs(line);
      if (line[0] == 'v')
	{
	  float x, y, z;
	  sscanf (line, "v %f %f %f", &x, &y, &z);
	  //printf("vertex xy: %f %f ", x, y);
	  if (minusx >= x)
	    minusx = x;
	  if (minusy >= y)
	    minusy = y;
	  if (minusz >= z)
	    minusz = z;
	  if (maxix <= x)
	    maxix = x;
	  if (maxiy <= y)
	    maxiy = y;
	  if (maxiz <= z)
	    maxiz = z;
	  if (minix >= x)
	    minix = (x > 0) ? x : 0;
	  if (miniy >= y)
	    miniy = (y > 0) ? y : 0;
	  if (miniz >= z)
	    miniz = (z > 0) ? z : 0;

	  vertexes[i][0] = x;
	  vertexes[i][1] = y;
	  vertexes[i][2] = z;

	  i++;
	}
      else if (line[0] == 'f')
	{
	  sscanf (line, "f %i %i %i", tray[j], tray[j] + 1, tray[j] + 2);
	  //printf("%i %i %i\n", tray[j][0], tray[j][1], tray[j][2]);
	  j++;
	}
    }
  size_v = i;
  size_tray = j;
  fclose (obj_file);


  float sf = 0.f;
  sf =
    (maxix - minusx >=
     maxiy - minusy) ? 1920.f / (maxix - minusx) : 1080.f / (maxiy - minusy);
  printf ("Type scale factor (recomended = %f): ", sf);
  scanf ("%f", &sf);

  float rx = 0, ry = 0, rz = 0;
  printf ("Type Rx rotation RADIANS( 0 = no operation ):\nRx: ");
  scanf ("%f", &rx);
  printf ("Type Ry rotation RADIANS( 0 = no operation ):\nRy: ");
  scanf ("%f", &ry);
  printf ("Type Rz rotation RADIANS( 0 = no operation ):\nRz: ");
  scanf ("%f", &rz);

  struct point
  {
    float x;
    float y;
    float z;
  } center =
    { (maxix + minusx + minix) / 2.f, (maxiy + minusy + miniy) / 2.f,
(maxiz + minusz + miniz) / 2.f };


  float tx = 0, ty = 0, tz = 0;
  printf
    ("Do you want to translate the object? format -> x y z: 0 0 0\nx y z: ");
  scanf ("%f %f %f", &tx, &ty, &tz);

  /*float f; */
  /*printf("Perspective proyection\nf="); */
  /*scanf("%f", &f); */

  // x rotation
  float Rx[4][4] = {
    {1, 0, 0, 0},
    {0, cos (rx), -sin (rx), 0},
    {0, sin (rx), cos (rx), 0},
    {0, 0, 0, 1}
  };
  // y rotation
  float Ry[4][4] = {
    {cos (ry), 0, -sin (ry), 0},
    {0, 1, 0, 0},
    {sin (ry), 0, cos (ry), 0},
    {0, 0, 0, 1}
  };
  // z rotation
  float Rz[4][4] = {
    {cos (rz), -sin (rz), 0, 0},
    {sin (rz), cos (rz), 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}
  };
  // translate to origin
  float To[4][4] = {
    {1, 0, 0, -center.x},
    {0, 1, 0, -center.y},
    {0, 0, 1, -center.z},
    {0, 0, 0, 1}
  };
  // translate to original pos
  float Tp[4][4] = {
    {1, 0, 0, center.x},
    {0, 1, 0, center.y},
    {0, 0, 1, center.z},
    {0, 0, 0, 1}
  };
  // scale factor 
  float S[4][4] = {
    {sf, 0, 0, 0},
    {0, sf, 0, 0},
    {0, 0, sf, 0},
    {0, 0, 0, 1}
  };

  float T[4][4] = {
    {1, 0, 0, tx},
    {0, 1, 0, ty},
    {0, 0, 1, tz},
    {0, 0, 0, 1}
  };

  float Tm[4][4] = {
    {1, 0, 0, -minusx - minix},
    {0, 1, 0, -minusy - miniy},
    {0, 0, 1, -minusz - miniz},
    {0, 0, 0, 1}
  };

  // Final matrix 
  float F[4][4] = {
    {1, 1, 1, 1},
    {1, 1, 1, 1},
    {1, 1, 1, 1},
    {1, 1, 1, 1},
  };

  // Auxiliary matrices for multiplication
  float aux1[4][4] = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1},
  }, aux2[4][4];

  // Perspective Proyection HAVE TROUBLES TROUBLESA TORBUELSFADFK
  /*float Pp[4][4] = { */
  /*{f,0,0,0}, */
  /*{0,f,0,0}, */
  /*{0,0,f,0}, */
  /*{0,0,1,0}, */
  /*}; */

  // Perform all matrix multiplications for the frists movements
  matrix_multiplication (Rx, To, aux1);	// move to origin and rotate on x
  matrix_multiplication (Ry, aux1, aux2);	// roration on y
  matrix_multiplication (Rz, aux2, aux1);	// rotation on z
  matrix_multiplication (Tp, aux1, aux2);	// move to original position
  matrix_multiplication (T, aux2, aux1);	// translate to somewhere
  matrix_multiplication (Tm, aux1, aux2);	// translate to positive coordinates
  // matrix_multiplication(Pp, aux2, aux1); // Perspective proyection TROUBLES
  matrix_multiplication (S, aux2, F);	// scale and save in final matrix 
  //matrix_multiplication(S, Tm, F);
  //matrix_multiplication(S, aux1, F);

  new_vertex (F);
  printf ("New Vertexes\n");

  minusx = 0, minusy = 0, minusz = 0, maxix = -FLT_MAX, maxiy =
    -FLT_MAX, maxiz = -FLT_MAX, miniy = FLT_MAX, minix = FLT_MAX, miniz =
    FLT_MAX;
  for (int v = 0; v < size_v; ++v)
    {
      if (maxix <= vertexes[v][0])
	maxix = vertexes[v][0];
      if (maxiy <= vertexes[v][1])
	maxiy = vertexes[v][1];
      if (maxiz <= vertexes[v][2])
	maxiz = vertexes[v][2];
    }

  // Pixel matrix for ppm file
  matrix = (int **) malloc (((int) maxiy + 1) * sizeof (int *));
  for (i = 0; i < (int) maxiy + 1; ++i)
    {
      matrix[i] = (int *) malloc (((int) maxix + 1) * sizeof (int));
    }

  Bresenham ();
  Save_ppm ("output.ppm");

  return 0;
}

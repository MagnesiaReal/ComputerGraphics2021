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


// Args for Drawings Methods with threads
struct args
{
  int ini, end;
};

struct xyz {
  float x;
  float y; 
  float z;
} show_all_planes;

struct pixel_array {
  int x,y;
  struct pixel_array *next;
};

struct edge {
  unsigned int v1;
  unsigned int v2;
} *edges;

struct face {
  unsigned int e1;
  unsigned int e2;
  unsigned int e3;
  struct pixel_array *pixel_array;
  unsigned int size_px;
} *faces;

char bufferus[128];
FILE *ppm_file;
double **vertexes, **vertexesd; // vertex array, vertex array copy// list of faces
int **matrix; // pixel 2D array for ppm file
int size_v = 0, size_faces =0, size_edges = 0; // size of vertex 
float minusx = 0, minusy = 0, minusz = 0, maxix = -FLT_MAX, maxiy =
  -FLT_MAX, maxiz = -FLT_MAX, miniy = FLT_MAX, minix = FLT_MAX, miniz =
  FLT_MAX;
pthread_t *threads;


void pixel_edge_list (int x, int y, int index) {
  struct pixel_array *ptr = (struct pixel_array *)malloc(sizeof (struct pixel_array));
  ptr->x = x;
  ptr->y = y;
  ptr->next = faces[index].pixel_array;
  faces[index].pixel_array = ptr;
}


// ##############################################################
void
Bresenham_draw_line (int x1, int y1, int x2, int y2, int i_edge)
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

  while (x != x2)
    {
      //printf("matrix[%i][%i]", x, y);
      //fflush(stdout);
      if ((*b >= 0 && *a >= 0) && (*b <= show_all_planes.y && *a <= show_all_planes.x))
	matrix[*b][*a] = 99;
      pixel_edge_list(*a, *b, i_edge);

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

void
Bresenham ()
{
  for (int i = 0; i < size_faces; ++i)
    {
      Bresenham_draw_line (vertexes[edges[faces[i].e1 - 1].v1 - 1][0],
			   vertexes[edges[faces[i].e1 - 1].v1 - 1][1],
			   vertexes[edges[faces[i].e1 - 1].v2 - 1][0],
			   vertexes[edges[faces[i].e1 - 1].v2 - 1][1], i);
      Bresenham_draw_line (vertexes[edges[faces[i].e2 - 1].v1 - 1][0],
			   vertexes[edges[faces[i].e2 - 1].v1 - 1][1],
			   vertexes[edges[faces[i].e2 - 1].v2 - 1][0],
			   vertexes[edges[faces[i].e2 - 1].v2 - 1][1], i);
      Bresenham_draw_line (vertexes[edges[faces[i].e3 - 1].v1 - 1][0],
			   vertexes[edges[faces[i].e3 - 1].v1 - 1][1],
			   vertexes[edges[faces[i].e3 - 1].v2 - 1][0],
			   vertexes[edges[faces[i].e3 - 1].v2 - 1][1], i);
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
      Bresenham_draw_line (vertexes[edges[faces[i].e1 - 1].v1 - 1][0],
			   vertexes[edges[faces[i].e1 - 1].v1 - 1][1],
			   vertexes[edges[faces[i].e1 - 1].v2 - 1][0],
			   vertexes[edges[faces[i].e1 - 1].v2 - 1][1], i);
      Bresenham_draw_line (vertexes[edges[faces[i].e2 - 1].v1 - 1][0],
			   vertexes[edges[faces[i].e2 - 1].v1 - 1][1],
			   vertexes[edges[faces[i].e2 - 1].v2 - 1][0],
			   vertexes[edges[faces[i].e2 - 1].v2 - 1][1], i);
      Bresenham_draw_line (vertexes[edges[faces[i].e3 - 1].v1 - 1][0],
			   vertexes[edges[faces[i].e3 - 1].v1 - 1][1],
			   vertexes[edges[faces[i].e3 - 1].v2 - 1][0],
			   vertexes[edges[faces[i].e3 - 1].v2 - 1][1], i);
      //printf("Figure %i complete\n", i);
      i++;
    }
  return 0;
}

// ##################################################3

void scan_line() {
  for (int fc = 0; fc < size_faces; ++fc) {
    for (struct pixel_array *px_itr = faces[fc].pixel_array; px_itr != NULL; px_itr = faces[fc].pixel_array->next) {
      
    }
  }
}


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
	      nv[i] += final_matrix[i][k] * vertexesd[v][k];
	    }
	}
      // save new vertex
      for (int i = 0; i < 4; ++i)
	{
	  vertexesd[v][i] = nv[i] / nv[3]; // Transform to no homogenious
	  vertexes[v][i] = nv[i] / nv[3];	// copy that
	  nv[i] = 0;
	}
    }
}

void Pp_vertex (float Pp_matrix[4][4]) {
  float nv[4] = { 0, 0, 0, 0 };
  for (int v = 0; v < size_v; ++v)
    {
      // vertex operation
      for (int i = 0; i < 4; ++i)
	{
	  for (int k = 0; k < 4; ++k)
	    {
	      nv[i] += Pp_matrix[i][k] * vertexes[v][k];
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
Save_ppm (char *filename)
{
  ppm_file = fopen (filename, "w");
  if (ppm_file == NULL)
    {
      perror ("No se pudo crear el archivo : ");
      exit (-1);
    }
  fprintf (ppm_file, "P3\n%i %i\n255", (int) (show_all_planes.x + 1), (int) (show_all_planes.y + 1));
  for (int i = show_all_planes.y; i >= 0; --i)
    {
      fputs ("\n", ppm_file);
      for (int j = 0; j <= show_all_planes.x; ++j)
	{
	  if (matrix[i][j] != 0)
	    fprintf (ppm_file, "%i %i %i ", matrix[i][j] + 61,
		     matrix[i][j] - 19, matrix[i][j] + 41);
	  else
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
      printf ("Use:\n./animation <vlf_file>");
      exit (EXIT_FAILURE);
    }
  obj_file = fopen (argv[1], "r");
  if (obj_file == NULL)
    {
      perror ("No se pudo abrir el archivo ");
      exit (-1);
    }
  // Building arrays of arrays for vertexes (x, y) and trays (v1, v2, v3)
  vertexes = (double **) malloc (Ver * sizeof (double *));
  for (int itr = 0; itr < Ver; itr++)
    {
      vertexes[itr] = (double *) malloc (4 * sizeof (double));
      vertexes[itr][3] = 1;
    }

  vertexesd = (double **) malloc (Ver * sizeof (double *));
  for (int itr = 0; itr < Ver; itr++)
    {
      vertexesd[itr] = (double *) malloc (4 * sizeof (double));
      vertexesd[itr][3] = 1;
    }

  faces = (struct face *) calloc (Tra, sizeof(struct face));
  for (int itr = 0; itr < Tra; ++itr) {
    faces[itr].pixel_array = NULL;
  }

  edges = (struct edge *)calloc(Tra*3, sizeof(struct edge ));
  

  char line[128];
  int i = 0, j = 0, k = 0, flag = -1;
  while (fgets (line, sizeof (line), obj_file) != NULL)
    {

      if (!strncmp(line, "vertex", 6)) {
	flag = 1;
	continue;
      } else if (!strncmp(line, "edges", 5)) {
        flag = 2;
	continue;
      } else if (!strncmp(line, "faces", 5)) {
        flag = 3;
	continue;
      }
      //fputs(line);
      if (flag == 1)
	{
	  double x, y, z;
	  sscanf (line, "%lf %lf %lf", &x, &y, &z);
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

	  vertexesd[i][0] = x;
	  vertexesd[i][1] = y;
	  vertexesd[i][2] = z;

	  i++;
	}
      else if (flag == 2) {
	sscanf(line, "%u %u", &edges[k].v1, &edges[k].v2);
	k++;
      }
      else if (flag == 3)
	{
	  sscanf (line, "%u %u %u", &faces[j].e1, &faces[j].e2 , &faces[j].e3);
	  //printf("%i %i %i\n", tray[j][0], tray[j][1], tray[j][2]);
	  j++;
	}
    }
  size_v = i;
  size_edges = k;
  size_faces = j;
  fclose (obj_file);

  /*for (int ci = 0; ci < size_v; ++ci) {*/
    /*printf("%lf %lf %lf\n", vertexesd[ci][0], vertexesd[ci][1], vertexesd[ci][2]);*/
  /*}*/

  /*for (int ki = 0; ki < size_edges; ++ki) {*/
    /*printf("%u %u\n", edges[ki].v1, edges[ki].v2);*/
  /*}*/
  
  /*for (int ki = 0; ki < size_faces; ++ki) {*/
    /*printf("%u %u %u\n", faces[ki].e1, faces[ki].e2, faces[ki].e3);*/
  /*}*/

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
    {1, 0, 0, 0 /*- minix*/ },
    {0, 1, 0, 0 /*- miniy*/ },
    {0, 0, 1, -minusz /*- miniz*/ },
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
      if (minusx >= vertexes[v][0])
	minusx = vertexes[v][0];
      if (minusy >= vertexes[v][1])
	minusy = vertexes[v][1];
      if (minusz >= vertexes[v][2])
	minusz = vertexes[v][2];
      if (maxix <= vertexes[v][0])
	maxix = vertexes[v][0];
      if (maxiy <= vertexes[v][1])
	maxiy = vertexes[v][1];
      if (maxiz <= vertexes[v][2])
	maxiz = vertexes[v][2];
      if (minix >= vertexes[v][0])
	minix = (vertexes[v][0] > 0) ? vertexes[v][0] : 0;
      if (miniy >= vertexes[v][1])
	miniy = (vertexes[v][1] > 0) ? vertexes[v][1] : 0;
      if (miniz >= vertexes[v][2])
	miniz = (vertexes[v][2] > 0) ? vertexes[v][2] : 0;
    }
  
  show_all_planes.x = (-minusx > maxix) ? maxix-2*minusx : 2*maxix;
  show_all_planes.y = (-minusy > maxiy) ? maxiy-2*minusy : 2*maxiy;

  printf ("maxX=%f maxY=%f\n", maxix, maxiy);
  printf ("Recalculate max and minimums done\n");

  // Automatic perspective projection
  float f = 0.8f*miniz;
  printf("%f\n",f);
  // Perspective Proyection HAVE TROUBLES TROUBLESA TORBUELSFADFK
  float Pp[4][4] = { 
  {f,0,0,0}, 
  {0,f,0,0}, 
  {0,0,f,0}, 
  {0,0,1,0}, 
  };

  // translate all planes
  float Ts[4][4] = {
    {1, 0, 0, (-minusx > maxix) ? -minusx : maxix},
    {0, 1, 0, (-minusy > maxiy) ? -minusy : maxiy},
    {0, 0, 1, 0},
    {0, 0, 0, 1}
  };

  float Fx[4][4] ={0};

  matrix_multiplication(Ts, Pp, Fx);
  Pp_vertex(Fx);

  // New center of object 
  center.x = (maxix + minusx + minix) / 2.f;
  center.y = (maxiy + minusy + miniy) / 2.f;
  center.z = (maxiz + minusz + miniz) / 2.f;

  // Pixel matrix for ppm file
  matrix = (int **) malloc (((int) show_all_planes.y + 1) * sizeof (int *));
  for (i = 0; i < (int) show_all_planes.y + 1; ++i)
    {
      matrix[i] = (int *) malloc (((int) show_all_planes.x + 1) * sizeof (int));
    }

  // translate to origin
  float Ton[4][4] = {
    {1, 0, 0, -center.x},
    {0, 1, 0, -center.y},
    {0, 0, 1, -center.z},
    {0, 0, 0, 1}
  };

  const float ryv = -0.123;
  float Ryv[4][4] = {
    {cos (ryv), 0, -sin (ryv), 0},
    {0, 1, 0, 0},
    {sin (ryv), 0, cos (ryv), 0},
    {0, 0, 0, 1}
  };

  // translate to original pos
  float Tpn[4][4] = {
    {1, 0, 0, center.x},
    {0, 1, 0, center.y},
    {0, 0, 1, center.z},
    {0, 0, 0, 1}
  };

  float Tv[4][4] = {
    {1, 0, 0, 0},
    {0, 1, 0, -4},
    {0, 0, 1, 0},
    {0, 0, 0, 1}
  };

  matrix_multiplication (Ryv, Ton, aux1);	// tranlate center object to origin and multiply by Y rotation
  matrix_multiplication (Tpn, aux1, aux2);	// translate object to original position
  matrix_multiplication (Tv, aux2, F);	// translate object four pixels down every frame
  printf ("Final matrix animation complete\n");

  int threadus = 3;		// number of threads
  threads = (pthread_t *) malloc (threadus * sizeof (pthread_t));	// allocate new trads

  for (int animation = 0; animation < 240; ++animation)
    {
      printf ("Init Bresenham Parallel Method\n");
      //Bresenham();    

      int tinc = size_faces / threadus, tend = 0, tini = 0;
      struct args argus[threadus];	// thread arguments 

      for (i = 0; i < threadus; ++i)
	{
	  tend += tinc;
	  argus[i].ini = tini;
	  if (i + 1 == threadus)
	    argus[i].end = size_faces;
	  else
	    argus[i].end = tend;
	  pthread_create (threads + i, NULL, Bresenham_thread,
			  (void *) &argus[i]);
	  tini = tend;
	}
      for (i = 0; i < threadus; ++i)
	{			// wait for all threads to finish
	  pthread_join (threads[i], NULL);
	}


      // file name 
      if (animation < 10)
	sprintf (bufferus, "00%i", animation);
      else if (animation < 100)
	sprintf (bufferus, "0%i", animation);
      else
	sprintf (bufferus, "%i", animation);
      strcat (bufferus, ".ppm");

      printf ("Saving to a File\n");
      Save_ppm (bufferus);

      // Recalculate vertexes for each frame
      new_vertex (F);
      Pp_vertex(Fx);

      // Reset the pixel matrix
      for (int i = 0; i <= show_all_planes.y; ++i)
	{
	  for (int j = 0; j <= show_all_planes.x; ++j)
	    {
	      matrix[i][j] = 0;
	    }
	}

    }

  return 0;
}

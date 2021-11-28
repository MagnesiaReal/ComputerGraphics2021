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

// Args for Drawings Methods with threads
struct args
{
  int ini, end;
};

struct xyz {
  int x, y, z;
} show_all_planes;
int **matrix = NULL; // pixel 2D array for ppm file

unsigned int size_edges;
struct edge {
  unsigned int v1;
  unsigned int v2;
} *edges = NULL;

unsigned int size_faces;
struct face {
  unsigned int e1;
  unsigned int e2;
  unsigned int e3;
  
  unsigned int px_size, fullfill_size;
  struct xyz *px_array, *fullfill_px;
} *faces = NULL;

int size_v = 0; // size of vertex 
double **vertexes = NULL, **vertexesd = NULL; // vertex array, vertex array copy// list of faces

float minusx = 0, minusy = 0, minusz = 0, maxix = -FLT_MAX, maxiy =
  -FLT_MAX, maxiz = -FLT_MAX, miniy = FLT_MAX, minix = FLT_MAX, miniz =
  FLT_MAX;

void save_pixel_per_face(const int x, const int y, const int face_num){
  faces[face_num].px_size++;
  struct xyz *array = (struct xyz *)realloc(faces[face_num].px_array, faces[face_num].px_size * sizeof(struct xyz));
  if(array != NULL) {
    faces[face_num].px_array = array;
    faces[face_num].px_array[faces[face_num].px_size-1].x = x;
    faces[face_num].px_array[faces[face_num].px_size-1].y = y;
  } else {
    perror("Error realloc in save_pixel_per_face ");
    exit(EXIT_FAILURE);
  }
}


// ##############################################################
void
Bresenham_draw_line (int x1, int y1, int x2, int y2, const int i_face)
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
      //printf("{x : %i, x2: %i}\n", x, x2);
      save_pixel_per_face(*a, *b, i_face);

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
  for (unsigned int i = 0; i < size_faces; ++i)
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

// ###################################################
// ############ FULL FILL ############################
int quick_sort_select(struct xyz array[], int left, int right) {
  int pivot = array[left].y;
  while (1) {
    while (array[left].y < pivot) {
      left++;
    }

    while (array[right].y > pivot) {
      right--;
    }

    if(left >= right) {
      return right;
    }

    struct xyz temp = array[left];
    array[left] = array[right];
    array[right] = temp;

    left++;
    right--;

  }
}


void quick_sort(struct xyz array[], int left, int right){
  if (left < right) {
    //printf("Sorting array[%i] = %i array[%i] = %i\n", left, array[left].y, right, array[right].y);
    int index_select = quick_sort_select(array, left, right);
    quick_sort(array, left, index_select);
    quick_sort(array, index_select + 1, right);
  }
}


struct xyz* scan_line(struct xyz* array, unsigned int *current_size, int min, int max, int y){
  int idx = *current_size;
  *current_size += max - min;
  
  if (max - min > 0)
  { 
    array = (struct xyz *)realloc(array, *current_size*sizeof(struct xyz));
    //printf("{current_size : %i, max : %i, min : %i, y : %i}\n", *current_size, max, min, y);
    
    if (array != NULL) { // No errors
      while(idx < *current_size){
        array[idx].x = ++min; 
        array[idx].y = y;
        idx++;
      }
    } else {
      perror("Error fullfill_px realloc ");
      exit(EXIT_FAILURE);
    }
  }

  return array;

}


void full_fill_faces() 
{  
  for (unsigned int f = 0; f < size_faces; f ++) {
    
    quick_sort(faces[f].px_array, 0, faces[f].px_size - 1);
    
    struct xyz *array = faces[f].px_array;
    unsigned int size = faces[f].px_size;
    unsigned int idx = 0;
    int max = INT_MIN, min = INT_MAX, yp = array[idx].y;

    while (idx < size){
      if(yp == array[idx].y){
        //printf("from min : %i max : %i\n", min, max);
        if (max < array[idx].x) max = array[idx].x;
        if (min > array[idx].x) min = array[idx].x;
        //printf("to min : %i max : %i\n", min, max);
        idx++;
      } else {
        //printf("scan line : {min : %i, max : %i}\n",min,max);
        faces[f].fullfill_px = scan_line(faces[f].fullfill_px, &faces[f].fullfill_size, min, max, array[idx].y);
        max = INT_MIN, min = INT_MAX;
        yp = array[idx].y;
      }

    }
  }
}

void raster() {
  printf("Rasterising\n");
  for(unsigned int f = 0; f < size_faces; f++) {

    struct xyz *array = faces[f].fullfill_px;
    for (unsigned int px = 0; px < faces[f].fullfill_size; px++) {
      if ((array[px].y >= 0 && array[px].x >= 0) && (array[px].y <= show_all_planes.y && array[px].x <= show_all_planes.x))
        matrix[array[px].y][array[px].x] = 130;
    }

    faces[f].fullfill_size = 0;
    free(faces[f].fullfill_px);
    faces[f].fullfill_px = NULL;
  }

  for (unsigned int f = 0; f < size_faces; ++f) {
    struct xyz *array = faces[f].px_array;
    for (unsigned int px = 0; px < faces[f].px_size; px++) {
      if ((array[px].y >= 0 && array[px].x >= 0) && (array[px].y <= show_all_planes.y && array[px].x <= show_all_planes.x))
        matrix[array[px].y][array[px].x] = 255;
    }

    faces[f].px_size = 0;
    free(faces[f].px_array);
    faces[f].px_array = NULL;
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
  FILE* ppm_file = fopen (filename, "w");
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
      fprintf (ppm_file, "%i %i %i ", matrix[i][j], matrix[i][j], matrix[i][j]);
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
    perror ("Can't open file ");
    exit (EXIT_FAILURE);
  }

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


      //New element
      vertexesd = (double **)realloc(vertexesd, ++i*sizeof(double*));
      if(vertexesd != NULL) {
        vertexesd[i-1] = (double *)calloc(4, sizeof(double));
        vertexesd[i-1][3] = 1;
        vertexesd[i-1][0] = x;
        vertexesd[i-1][1] = y;
        vertexesd[i-1][2] = z;
      } else {
        perror("Error reallocating vertexesd");
        exit(EXIT_FAILURE);
      }
    }
    else if (flag == 2) {
      edges = (struct edge *)realloc(edges, ++k*sizeof(struct edge));
      if (edges != NULL) {
        sscanf(line, "%u %u", &edges[k-1].v1, &edges[k-1].v2);   
      } else {
        perror("Error reallocating edges ");
        exit(EXIT_FAILURE);
      }
    }
    else if (flag == 3)
    {
      faces = (struct face *)realloc(faces, ++j*sizeof(struct face));
      if (faces != NULL) {
        sscanf (line, "%u %u %u", &faces[j-1].e1, &faces[j-1].e2 , &faces[j-1].e3);
        faces[j-1].px_array = NULL; // 0
        faces[j-1].px_size = 0;
        faces[j-1].fullfill_px = NULL;
        faces[j-1].fullfill_size = 0;
      } else {
        perror("Error reallocating faces ");
        exit(EXIT_FAILURE);
      }
    }
  }
  size_v = i;
  size_edges = k;
  size_faces = j;

  fclose (obj_file);

  vertexes = (double **)malloc(i*sizeof(double*));
  if(vertexes != NULL) {
    for (int ir = 0; ir < i; ++ir) {
      vertexes[ir] = (double *)calloc(4, sizeof(double)); 
    }
  } else {
    perror("Error reallocating vertexes");
    exit(EXIT_FAILURE);
  }

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

  printf ("maxX=%i maxY=%i\n", show_all_planes.x, show_all_planes.y);
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
  pthread_t *threads = (pthread_t *) malloc (threadus * sizeof (pthread_t));	// allocate new trads

  for (int animation = 0; animation < 240; ++animation)
  {
    printf ("%i Init Bresenham Parallel Method\n", animation);
    // BRESENHAM /////////////////////////////////////////
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
    ////////////////////////////////////////////////////////
    // FULLFILL PXs ////////////////////////////////////////
    full_fill_faces();
    raster();

    // file name 
    char bufferus[128];
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

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistdio.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <pthread.h>

#define Ver 2001000
#define Tra 4000000

// params for Threads
struct args
{
  int ini, end;
};


FILE *ppm_file;
int **vertexes;
int **tray;
int **matrix;
float fact_scale = 10.f;
int minusx = 0, minusy = 0, maxix = INT_MIN, maxiy = INT_MIN, miniy =
  INT_MAX, minix = INT_MAX;
pthread_t *threads;
// ############################################################################
// ############################################################################
// ############################################################################
void
Naive_draw_line (int x1, int y1, int x2, int y2)
{
  float m = ((float) y2 - y1) / ((float) x2 - x1);
  //printf("x1=%i y1=%i x2=%i y2=%i m=%f\n", x1, y1, x2, y2, m);

  // if -nan
  if (x2 == x1)
    {
      int y, ymax;

      if (y1 <= y2)
	{
	  y = y1;
	  ymax = y2;
	}
      else
	{
	  y = y2;
	  ymax = y1;
	}

      while (y <= ymax)
	{
	  //printf("matrix[%i][%i] = 150 to equals X\n", maxiy - y, x1 - minusx - minix);
	  matrix[maxiy - y][x1 - minusx - minix] = 150;
	  y++;
	};
      return;
    }
  int x, xmax, y;
  if (x1 < x2)
    {
      x = x1;
      xmax = x2;
    }
  else
    {
      x = x2;
      xmax = x1;
    }
  while (x <= xmax)
    {
      y = m * (x - x1) + y1;
      //if (maxiy - y > miniy) printf("matrix[%i][%i] = 150\n", maxiy - y, x-minusx - minix);
      //if (maxiy - y > miniy) printf("x1=%i y1=%i x2=%i y2=%i m=%f\n", x1, y1, x2, y2, m);
      matrix[maxiy - y][x - minusx - minix] = 150;	// nivel de color
      x++;
    }
}

// Funcion de la forma y = m(x) + b    x1 < x < x2 || x2 < x < x1 
void
Naive ()
{
  for (int i = 0; i < Tra; ++i)
    {
      if (tray[i][0] == 0)
	break;
      Naive_draw_line (vertexes[tray[i][0] - 1][0],
		       vertexes[tray[i][0] - 1][1],
		       vertexes[tray[i][1] - 1][0],
		       vertexes[tray[i][1] - 1][1]);
      Naive_draw_line (vertexes[tray[i][1] - 1][0],
		       vertexes[tray[i][1] - 1][1],
		       vertexes[tray[i][2] - 1][0],
		       vertexes[tray[i][2] - 1][1]);
      Naive_draw_line (vertexes[tray[i][2] - 1][0],
		       vertexes[tray[i][2] - 1][1],
		       vertexes[tray[i][0] - 1][0],
		       vertexes[tray[i][0] - 1][1]);
      //printf("Figure %i complete\n", i);
    }
}

void *
Naive_thread (void *argsp)
{
  struct args *arguments = (struct args *) argsp;
  int i = arguments->ini, end = arguments->end;
  //printf("ini=%i end=%i\n", i, end);
  while (i < end)
    {
      if (tray[i][0] == 0)
	break;
      Naive_draw_line (vertexes[tray[i][0] - 1][0],
		       vertexes[tray[i][0] - 1][1],
		       vertexes[tray[i][1] - 1][0],
		       vertexes[tray[i][1] - 1][1]);
      Naive_draw_line (vertexes[tray[i][1] - 1][0],
		       vertexes[tray[i][1] - 1][1],
		       vertexes[tray[i][2] - 1][0],
		       vertexes[tray[i][2] - 1][1]);
      Naive_draw_line (vertexes[tray[i][2] - 1][0],
		       vertexes[tray[i][2] - 1][1],
		       vertexes[tray[i][0] - 1][0],
		       vertexes[tray[i][0] - 1][1]);
      //printf("Figure %i complete\n", i);
      i++;
    }
  return 0;
}

// ###########################################################################
// ############################################################################
// ############################################################################
void
DDA_draw_line (int x1, int y1, int x2, int y2)
{
  float m = ((float) y2 - y1) / ((float) x2 - x1);
  //printf("x1=%i y1=%i x2=%i y2=%i m=%f\n", x1, y1, x2, y2, m);

  // if -nan nothing change
  if (x2 == x1)
    {
      int y, ymax;
      if (y1 <= y2)
	{
	  y = y1;
	  ymax = y2;
	}
      else
	{
	  y = y2;
	  ymax = y1;
	}
      while (y <= ymax)
	{
	  if (maxiy - y > maxiy)
	    printf ("matrix[%i][%i] = 150\n", maxiy - y, x1 - minusx - minix);
	  matrix[maxiy - y][x1 - minusx - minix] = 150;
	  y++;
	};
      return;
    }
  int x, xmax, y;
  float yf;
  if (x1 < x2)
    {
      x = x1;
      xmax = x2;
    }
  else
    {
      x = x2;
      xmax = x1;
    }
  yf = m * (x - x1) + y1;
  y = yf;
  while (x <= xmax)
    {
      //printf("matrix[%i][%i] ", maxiy - y, x-minusx - minix);
      //fflush(stdout); // for print all statments before crash
      matrix[(maxiy - y > maxiy - minusy - miniy) ? maxiy - minusy - miniy : maxiy - y][x - minusx - minix] = 150;	// level color 150
      x++;

      yf = yf + m;		// young y = old y + m;
      y = yf;
    }
}

void
DDA ()
{
  for (int i = 0; i < Tra; ++i)
    {
      if (tray[i][0] == 0)
	break;
      DDA_draw_line (vertexes[tray[i][0] - 1][0], vertexes[tray[i][0] - 1][1],
		     vertexes[tray[i][1] - 1][0],
		     vertexes[tray[i][1] - 1][1]);
      DDA_draw_line (vertexes[tray[i][1] - 1][0], vertexes[tray[i][1] - 1][1],
		     vertexes[tray[i][2] - 1][0],
		     vertexes[tray[i][2] - 1][1]);
      DDA_draw_line (vertexes[tray[i][2] - 1][0], vertexes[tray[i][2] - 1][1],
		     vertexes[tray[i][0] - 1][0],
		     vertexes[tray[i][0] - 1][1]);
      //printf("Figure %i complete\n", i);
    }
}

void *
DDA_thread (void *argsp)
{
  struct args *arguments = (struct args *) argsp;
  int i = arguments->ini, end = arguments->end;
  //printf("ini=%i end=%i\n", i, end);
  while (i < end)
    {
      if (tray[i][0] == 0)
	break;
      DDA_draw_line (vertexes[tray[i][0] - 1][0], vertexes[tray[i][0] - 1][1],
		     vertexes[tray[i][1] - 1][0],
		     vertexes[tray[i][1] - 1][1]);
      DDA_draw_line (vertexes[tray[i][1] - 1][0], vertexes[tray[i][1] - 1][1],
		     vertexes[tray[i][2] - 1][0],
		     vertexes[tray[i][2] - 1][1]);
      DDA_draw_line (vertexes[tray[i][2] - 1][0], vertexes[tray[i][2] - 1][1],
		     vertexes[tray[i][0] - 1][0],
		     vertexes[tray[i][0] - 1][1]);
      //printf("Figure %i complete\n", i);
      i++;
    }
  return 0;
}

// ############################################################################
// ############################################################################
// ############################################################################

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
      matrix[maxiy - y][x - minusx - minix] = 150;

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

// ############################################################################
// ############################################################################
// ############################################################################


void
Save_ppm ()
{
  ppm_file = fopen ("model.ppm", "w");
  if (ppm_file == NULL)
    {
      perror ("No se pudo crear el archivo : ");
      exit (-1);
    }
  fprintf (ppm_file, "P3\n%i %i\n255", maxix - minusx - minix + 1,
	   maxiy - minusy - miniy + 1);
  for (int i = 0; i < maxiy - minusy - miniy + 1; ++i)
    {
      fputs ("\n", ppm_file);
      for (int j = 0; j < maxix - minusx - minix + 1; ++j)
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
  // Reading File
  FILE *obj_file;
  if (argc < 2)
    {
      printf ("Use:\n./objtoppm <objfile>");
      exit (EXIT_FAILURE);
    }
  obj_file = fopen (argv[1], "r");

  if (obj_file == NULL)
    {
      perror ("No se pudo abrir el archivo : ");
      exit (-1);
    }
  // Building arrays of arrays for vertexes (x, y) and trays (v1, v2, v3)
  vertexes = (int **) malloc (Ver * sizeof (int *));
  for (int itr = 0; itr < Ver; itr++)
    {
      vertexes[itr] = (int *) malloc (2 * sizeof (int));
    }

  tray = (int **) malloc (Tra * sizeof (int *));
  for (int itr = 0; itr < Tra; itr++)
    {
      tray[itr] = (int *) malloc (3 * sizeof (int));
    }

  // Scale could be 
  float min_sx = 1080.f / (maxix - minusx - minix),
    min_sy = 1920.f / (maxiy - minusy - miniy);

  if (min_sx <= min_sy)
    fact_scale = min_sx;
  else
    fact_scale = min_sy;

  printf ("Type scale: ");
  scanf ("%f", &fact_scale);


  char line[128];
  int i = 0, j = 0;
  while (fgets (line, sizeof (line), obj_file) != NULL)
    {
      //fputs(line);
      if (line[0] == 'v')
	{
	  float x, y;
	  sscanf (line, "v %f %f ", &x, &y);
	  //printf("vertex xy: %f %f ", x, y);
	  // Observamos el valor mas negativo que existe y lo guardamos
	  int xt = x * fact_scale;
	  int yt = y * fact_scale;
	  if (minusx >= xt)
	    minusx = xt;
	  if (minusy >= yt)
	    minusy = yt;
	  if (maxix <= xt)
	    maxix = xt;
	  if (maxiy <= yt)
	    maxiy = yt;
	  if (minix >= xt)
	    minix = (xt > 0) ? xt : 0;
	  if (miniy >= yt)
	    miniy = (yt > 0) ? yt : 0;
	  vertexes[i][0] = xt;
	  vertexes[i][1] = yt;
	  i++;
	}
      else if (line[0] == 'f')
	{
	  sscanf (line, "f %i %i %i", tray[j], tray[j] + 1, tray[j] + 2);
	  //printf("%i %i %i\n", tray[j][0], tray[j][1], tray[j][2]);
	  j++;
	}
    }
  fclose (obj_file);


  printf
    ("min positive value y (miniy) = %i\nmin positive value x (minix) = %i\n",
     miniy, minix);
  printf
    ("max negative value y (minusy) = %i\nmax negative value x (minusx) = %i\n",
     minusy, minusx);
  printf ("max value x (maxix) = %i\nmin value x (minusx) = %i\n", maxix,
	  minusx);
  printf ("max value y (maxiy) = %i\nmin value y (minusy) = %i\n", maxiy,
	  minusy);
  printf ("min wieght (maxix - minusx - minix) = %i\n",
	  maxix - minusx - minix);
  printf ("min height (maxiy - minusy - miniy) = %i\n",
	  maxiy - minusy - miniy);

  matrix = (int **) malloc ((maxiy - minusy + 1 - miniy) * sizeof (int *));
  for (i = 0; i < maxiy - minusy + 1 - miniy; ++i)
    {
      matrix[i] =
	(int *) malloc ((maxix - minusx + 1 - minix) * sizeof (int));
    }

  //DDA_draw_line(-1988, 17008, -2007, 17009);

  int option = 1, threadus, tini = 0, tend = 0;

  printf ("\nChoise Method to draw XY proyection:\n");
  printf
    ("1 Native\n2 DDA\n3 Bresenham\n4 Parallel Naive\n5 Parallel DDA\n6 Parallel Bresenham\noption: ");
  fflush (stdin);
  scanf ("%d", &option);
  fflush (stdin);

  if (option > 3 && option < 7)
    {
      printf ("How many threads do you want?\n");
      printf ("threads: ");
      scanf ("%i", &threadus);
      threads = (pthread_t *) malloc (threadus * sizeof (pthread_t));
    }

  clock_t startus = clock (), endus;
  switch (option)
    {
    case 1:
      printf ("Init Naive Method\n");
      Naive ();
      break;
    case 2:
      printf ("Init DDA Method\n");
      DDA ();
      break;
    case 3:
      printf ("Init Bresenham Method\n");
      Bresenham ();
      break;
    }

  if (option == 4)
    {
      printf ("Init Naive Parallel Method\n");
      int tinc = Tra / threadus;
      struct args argus[threadus];
      for (i = 0; i < threadus; ++i)
	{
	  tend += tinc;
	  argus[i].ini = tini;
	  if (i + 1 == threadus)
	    argus[i].end = Tra;
	  else
	    argus[i].end = tend;
	  pthread_create (threads + i, NULL, Naive_thread,
			  (void *) &argus[i]);
	  tini = tend;
	}
      for (i = 0; i < threadus; ++i)
	{
	  pthread_join (threads[i], NULL);
	}
    }
  else if (option == 5)
    {
      printf ("Init DDA Parallel Method\n");
      int tinc = Tra / threadus;
      struct args argus[threadus];
      for (i = 0; i < threadus; ++i)
	{
	  tend += tinc;
	  argus[i].ini = tini;
	  if (i + 1 == threadus)
	    argus[i].end = Tra;
	  else
	    argus[i].end = tend;
	  pthread_create (threads + i, NULL, DDA_thread, (void *) &argus[i]);
	  tini = tend;
	}
      for (i = 0; i < threadus; ++i)
	{
	  pthread_join (threads[i], NULL);
	}
    }
  else if (option == 6)
    {
      printf ("Init Bresenham Parallel Method\n");
      int tinc = Tra / threadus;
      struct args argus[threadus];
      for (i = 0; i < threadus; ++i)
	{
	  tend += tinc;
	  argus[i].ini = tini;
	  if (i + 1 == threadus)
	    argus[i].end = Tra;
	  else
	    argus[i].end = tend;
	  pthread_create (threads + i, NULL, Bresenham_thread,
			  (void *) &argus[i]);
	  tini = tend;
	}
      for (i = 0; i < threadus; ++i)
	{
	  pthread_join (threads[i], NULL);
	}
    }

  endus = clock ();
  double time_method = (double) (endus - startus) / CLOCKS_PER_SEC;
  printf ("Method complete in %f ms\n", time_method);
  printf ("Saving PPM file\n");
  Save_ppm ();
}

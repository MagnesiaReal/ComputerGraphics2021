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
#include <stdbool.h>

#define Ver 2001000
#define Tra 4000000

struct v_xyz {
  double x, y, z, w;
} *vertexesd = NULL;

struct edge {
  unsigned int v1, v2;
} *edges;

struct face {
  unsigned int e1, e2, e3;
} *faces; 

int main(int argc, char *argv[])
{

  FILE *obj_file;
  if (argc < 2)
  {
    printf ("Use:\n./objtovlf <obj_file>");
    exit (EXIT_FAILURE);
  }
  
  obj_file = fopen (argv[1], "r");
  if (obj_file == NULL)
  {
    perror ("No se pudo abrir el archivo ");
    exit (-1);
  }
 
  /////////////////////////////////////////////////////////////////
  char line[128];
  int i = 0, j = 0, e = 0;
  while (fgets (line, sizeof (line), obj_file) != NULL)
  {     
    //fputs(line);
    if (line[0] == 'v')
    {
      vertexesd = (struct v_xyz*)realloc(vertexesd, (i+1)*sizeof(struct v_xyz));

      if (vertexesd != NULL) {
        double x, y, z;
        sscanf (line, "v %lf %lf %lf", &vertexesd[i].x, &vertexesd[i].y, &vertexesd[i].z);
        vertexesd[i].w = 1;
        i++;
      } else {
        perror("Error in vertexes realloc ");
        exit(EXIT_FAILURE);
      }
    }
    else if (line[0] == 'f')
    {
      edges = (struct edge*)realloc(edges, (e+3)*sizeof(struct edge));
      faces = (struct face*)realloc(faces, (j+1)*sizeof(struct face));

      if (edges != NULL && faces != NULL) {
        unsigned int vexn[3];
        sscanf (line, "f %i %i %i", vexn, vexn+1, vexn+2);
        edges[e].v1 = vexn[0];
        edges[e].v2 = vexn[1];
        faces[j].e1 = e + 1;
        e++;

        edges[e].v1 = vexn[1];
        edges[e].v2 = vexn[2];
        faces[j].e2 = e + 1;
        e++;

        edges[e].v1 = vexn[2];
        edges[e].v2 = vexn[0];
        faces[j].e3 = e + 1;
        e++;
        j++;
      } else {
        perror("Error in edges or faces ");
        exit(EXIT_FAILURE);
      }    
    }
  }
  fclose (obj_file);
  /////////////////////////////////////////////////////////////////
  char name[128];
  sprintf(name, "%s.vlf",argv[1]);

  FILE* vlf_file = fopen(name, "w");

  fprintf(vlf_file, "vertex %i\n", i);
  for (int v = 0; v < i; ++v) {
    fprintf(vlf_file, "%lf %lf %lf\n", vertexesd[v].x, vertexesd[v].y, vertexesd[v].z);
  }
  // getting faces and edges
  fprintf(vlf_file, "edges %i\n", e);
  for (int er = 0; er < e; ++er) {
    fprintf(vlf_file, "%u %u\n", edges[er].v1, edges[er].v2);
  }

  fprintf(vlf_file, "faces %i\n", j);
  for (int fa = 0; fa < j; ++fa) {
    fprintf(vlf_file, "%u %u %u\n", faces[fa].e1, faces[fa].e2, faces[fa].e3);
  }

  fclose(vlf_file);
  return 0;
}

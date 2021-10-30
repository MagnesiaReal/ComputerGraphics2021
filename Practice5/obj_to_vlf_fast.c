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

double **vertexesd;
unsigned int **edges, **faces, **facesinc;
int i = 0, j = 0, e = 0, checkpoint = INT_MAX;

void verif(int *ir, const int v1, const int v2, int start) {
  for (int er = start + 1; er < *ir; ++er) {
    if ((edges[er][0] == v1 && edges[er][1] == v2) || (edges[er][0] == v2 && edges[er][1] == v1)) {
      printf("find repeat: %u %u\n", edges[er][0], edges[er][1]);
      for (int k = er; k < *ir-1; ++k) {
	edges[k][0] = edges[k+1][0];
	edges[k][1] = edges[k+1][1];
      }

      for (int k = start + 1; k < e; ++k) {
        if(faces[(k)/3][(k)%3] > er)
	    faces[(k)/3][(k)%3]--;
      }

      /*printf("er = %i\n", er);*/
      /*int auxinc = facesinc[er/3][er%3];*/
      /*printf("preaux = %i\n", auxinc);*/

      /*auxinc = facesinc[(er+auxinc)/3][(er+auxinc)%3];*/
      /*printf("faceinc = %i\n", auxinc);*/
      
      /*while (start + 1 > faces[(er+auxinc)/3][(er+auxinc)%3]) {*/
	/*auxinc = facesinc[(er+auxinc)/3][(er+auxinc)%3];*/
	/*printf("new auxinc = %i\n", auxinc);*/
      /*}	*/

      /*if(checkpoint > er) {*/
	/*for (int kr = er; kr < e; ++kr) {*/
	  /*if(faces[(kr)/3][(kr)%3] > er)*/
	    /*faces[(kr)/3][(kr)%3]--;*/
	  /*facesinc[kr/3][kr%3] += 1;*/
	/*}*/
	/*faces[(er)/3][(er)%3] = start + 1;*/
	/*printf("er = %i\n",(er));*/
	/*checkpoint = er;*/
      /*} else {*/
	/*for (int kr = er + auxinc; kr < e; ++kr) {*/
	  /*if(faces[(kr)/3][(kr)%3] > er)*/
	    /*faces[(kr)/3][(kr)%3]--;*/
	  /*facesinc[(kr)/3][(kr)%3] += 1;*/
	  /*printf("faceinc[%i] = %i\n", kr, facesinc[kr/3][kr%3]);*/
	/*}*/
	/*faces[(er+auxinc)/3][(er+auxinc)%3] = start + 1;*/
	/*printf("er+aux = %i\n",(er+auxinc));*/
      /*}*/

      /*for (int jar = 0; jar < j; ++jar) {*/
	/*printf("%u %u %u\n",faces[jar][0], faces[jar][1], faces[jar][2]);*/
      /*}*/

      *ir = *ir - 1;
      er--;
    }
  }
}


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
  // 2D Arrays ////////////////////////////////////////////////////
  vertexesd = (double **) malloc (Ver * sizeof (double *));
  for (int itr = 0; itr < Ver; itr++)
    {
      vertexesd[itr] = (double *) malloc (4 * sizeof (double));
      vertexesd[itr][3] = 1;
    }

  faces = (unsigned int **) malloc (Tra * sizeof (unsigned int *));
  for (int itr = 0; itr < Tra; itr++)
    {
      faces[itr] = (unsigned int *) malloc (3 * sizeof (unsigned int));
    }
  
  edges = (unsigned int **)malloc(Tra*3*sizeof(unsigned int *));
  for (int it = 0; it < Tra*3; it++){
    edges[it] = (unsigned int *)malloc(2*sizeof(unsigned int));
  }
  /////////////////////////////////////////////////////////////////
  char line[128];
  
  while (fgets (line, sizeof (line), obj_file) != NULL)
    {
      //fputs(line);
      if (line[0] == 'v')
	{
	  double x, y, z;
	  sscanf (line, "v %lf %lf %lf", &x, &y, &z);
	  //printf("vertex xy: %f %f ", x, y);
	  vertexesd[i][0] = x;
	  vertexesd[i][1] = y;
	  vertexesd[i][2] = z;

	  i++;
	}
      else if (line[0] == 'f')
	{
	  unsigned int vexn[3];
	  sscanf (line, "f %i %i %i", vexn, vexn+1, vexn+2);
	  edges[e][0] = vexn[0];
	  edges[e][1] = vexn[1];
	  faces[j][0] = e + 1;
	  e++;

	  edges[e][0] = vexn[1];
	  edges[e][1] = vexn[2];
	  faces[j][1] = e + 1;
	  e++;

	  edges[e][0] = vexn[2];
	  edges[e][1] = vexn[0];
	  faces[j][2] = e + 1;
	  e++;

	  //printf("face number %i\n", j);
	  j++;
	}
    }
  fclose (obj_file);
  int sizer = e;
  

  // Change repeat edges for first find
  printf("Searching repeat lines\n");
  for (int ie = 0; ie < e; ++ie) {
    printf("itr %i\n", ie);
    unsigned int v1 = edges[ie][0], v2 = edges[ie][1];
    for (int je = ie + 1; je < e; ++je) {
      if ((edges[je][0] == v1 && edges[je][1] == v2) || (edges[je][0] == v2 && edges[je][1] == v1)) {
	faces[je/3][je%3] = ie + 1;
      }
    }
  }


  // Trying to remove duplicated lines;
  printf("Removing repeat lines\n");
  for (int iar = 0; iar < sizer; ++iar) {
    verif(&sizer, edges[iar][0], edges[iar][1], iar);
  }

  
  char name[128];
  sprintf(name, "%s.vlf",argv[1]);
  FILE* vlf_file = fopen(name, "w");
  
  fprintf(vlf_file, "vertex %i\n", i);
  for (int v = 0; v < i; ++v) {
    fprintf(vlf_file, "%lf %lf %lf\n", vertexesd[v][0], vertexesd[v][1], vertexesd[v][2]);
    /*fwrite(&vertexesd[v][0], 1, sizeof(vertexesd[v][0]), vlf_file);*/
    /*fwrite(&vertexesd[v][1], 1, sizeof(vertexesd[v][1]), vlf_file);*/
    /*fwrite(&vertexesd[v][2], 1, sizeof(vertexesd[v][2]), vlf_file);*/
  }
  // getting faces and edges
  fprintf(vlf_file, "edges %i\n", sizer);
  for (int er = 0; er < sizer; ++er) {
    fprintf(vlf_file, "%u %u\n", edges[er][0], edges[er][1]);
    /*fwrite(&edges[er][0], 1, sizeof(edges[er][0]), vlf_file);*/
    /*fwrite(&edges[er][1], 1, sizeof(edges[er][1]), vlf_file);*/
  }
  
  fprintf(vlf_file, "faces %i\n", j);
  for (int fa = 0; fa < j; ++fa) {
    fprintf(vlf_file, "%u %u %u\n", faces[fa][0], faces[fa][1], faces[fa][2]);
    /*fwrite(&faces[fa][0], 1, sizeof(faces[fa][0]), vlf_file);*/
    /*fwrite(&faces[fa][1], 1, sizeof(faces[fa][1]), vlf_file);*/
    /*fwrite(&faces[fa][2], 1, sizeof(faces[fa][2]), vlf_file);*/
  }

  fclose(vlf_file);
  return 0;
}

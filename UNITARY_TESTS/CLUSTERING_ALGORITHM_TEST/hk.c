/* Tobin Fricke's implementation of the Hoshen-Kopelman algorithm for
   cluster labeling (library version, no main function here).

   Copyright (c) September 9, 2000, by Tobin Fricke <tobin@splorg.org>
   Distributed under the terms of the GNU Public License.

   Modified 2002-03-09 Tobin Fricke
   Modified substantially 2004-04-21 by Tobin Fricke

   This program is written in the 1999 standard of the C language (C99).  Older C
   compilers will refuse to compile it.   You can use a C++ compiler, a C99 compiler,
   or you can modify this code to comply with a previous version of the C standard.
   The GCC compiler supports C99 as of version 3.0.  Compile the library with:

   gcc-3.0 -Wall -std=c99 -c hk.c

   http://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "hk.h"

/* Implementation of Union-Find Algorithm (as a library of functions)
   There is no main() function here: see Hoshen-Kopelman.c file for a full test 
   main program.
*/

/* The 'labels' array has the meaning that labels[x] is an alias for the label x; by
   following this chain until x == labels[x], you can find the canonical name of an
   equivalence class.  The labels start at one; labels[0] is a special value indicating
   the highest label already used. */

int *labels;
int  n_labels = 0;     /* length of the labels array */

/*  uf_find returns the canonical label for the equivalence class containing x */
int uf_find(int x) {
  int y = x;
  while (labels[y] != y)
    y = labels[y];
  
  while (labels[x] != x) {
    int z = labels[x];
    labels[x] = y;
    x = z;
  }
  return y;
}

/*  uf_union joins two equivalence classes and returns the canonical label of the resulting class. */
int uf_union(int x, int y) {
  return labels[uf_find(x)] = uf_find(y);
}

/*  uf_make_set creates a new equivalence class and returns its label */
int uf_make_set(void) {
  labels[0] ++;
  assert(labels[0] < n_labels);
  labels[labels[0]] = labels[0];
  return labels[0];
}

/*  uf_intitialize sets up the data structures needed by the union-find implementation. */
void uf_initialize(int max_labels) {
  n_labels = max_labels;
  labels = calloc(sizeof(int), n_labels);
  labels[0] = 0;
}

/*  uf_done frees the memory used by the union-find data structures */
void uf_done(void) {
  n_labels = 0;
  free(labels);
  labels = 0;
}
/* End Union-Find implementation */

#define max(a,b) (a>b?a:b)
#define min(a,b) (a>b?b:a)

/* print_matrix prints out a matrix that is set up in the "pointer to pointers" scheme
   (aka, an array of arrays); this is incompatible with C's usual representation of 2D
   arrays, but allows for 2D arrays with dimensions determined at run-time */
void print_matrix(int **matrix, int m, int n) {
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++)
      printf("%3d ",matrix[i][j]);
    printf("\n");
  }
}

/* Label the clusters in "matrix".  Return the total number of clusters found. */
int hoshen_kopelman(int **matrix, int m, int n) {
  
  uf_initialize(m * n / 2);
  
  /* scan the matrix */
  
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      if (matrix[i][j]) {                        // if occupied ...

	      int up = (i==0 ? 0 : matrix[i-1][j]);    //  look up  
	      int left = (j==0 ? 0 : matrix[i][j-1]);  //  look left
	
	      switch (!!up + !!left) {
	  
	      case 0:
	        matrix[i][j] = uf_make_set();      // a new cluster
	      break;
	  
	      case 1:                              // part of an existing cluster
	        matrix[i][j] = max(up,left);       // whichever is nonzero is labelled
	      break;
	  
	      case 2:                              // this site binds two clusters
	        matrix[i][j] = uf_union(up, left);
	      break;
	    }
    }
  
  /* apply the relabeling to the matrix */

  /* This is a little bit sneaky.. we create a mapping from the canonical labels
     determined by union/find into a new set of canonical labels, which are 
     guaranteed to be sequential. */
  
  int *new_labels = calloc(sizeof(int), n_labels); // allocate array, initialized to zero
  
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      if (matrix[i][j]) {
	      int x = uf_find(matrix[i][j]);
	        if (new_labels[x] == 0) {
	          new_labels[0]++;
	          new_labels[x] = new_labels[0];
	        }
	      matrix[i][j] = new_labels[x];
      }
 
  int total_clusters = new_labels[0];

  free(new_labels);
  uf_done();

  return total_clusters;
}

/* This procedure checks to see that any occupied neighbors of an occupied site
   have the same label. */

void check_labelling(int **matrix, int m, int n) 
{
  int N,S,E,W;
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      if (matrix[i][j]) {
	      N = ( i==0 ? 0 : matrix[i-1][j] );
	      S = ( i==m-1 ? 0 : matrix[i+1][j] );
	      E = ( j==n-1 ? 0 : matrix[i][j+1] );
	      W = ( j==0 ? 0 : matrix[i][j-1] );
	
	      assert( N==0 || matrix[i][j]==N );
	      assert( S==0 || matrix[i][j]==S );
	      assert( E==0 || matrix[i][j]==E );
	      assert( W==0 || matrix[i][j]==W );
      }
}

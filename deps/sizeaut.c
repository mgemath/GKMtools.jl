/*
  This code is adapted from the program nautyex1.c in the nauty distribution.
  The input is an adjacency matrix of the graph, the order of the matrix is n.
  The output is the size of the automorphism group.
  The graph is assumed to have at most MAXN vertices.
*/

#define MAXN 1000 //maximum number of vertices
#include <math.h>
#include <nauty.h>


unsigned long int sizeaut(int *mat, long int n)
{
  long unsigned int g[MAXN*MAXM];
  int lab[MAXN],ptn[MAXN],orbits[MAXN];
  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;
  options.writeautoms = FALSE;
  int m = SETWORDSNEEDED(n);
  EMPTYGRAPH(g,m,n);
  for (int v = 0; v < (n-1); ++v) 
  {
    for (int u = (v+1); u < n; ++u) 
    {
      if (mat[v*n + u])
      {
        ADDONEEDGE(g,v,u,m);
      }
    }
  }
  
  densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);
  return (stats.grpsize1 * (unsigned long) pow(10, stats.grpsize2));
}
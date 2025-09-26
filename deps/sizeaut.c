#include <math.h>
#include <nauty.h>
//g is a matrix n*n, n must be < 1000
// int sizeaut(long unsigned int *g, int n)
unsigned long int sizeaut(int *mat, long int n)
{
  long unsigned int g[n*n];
  int lab[1000],ptn[1000],orbits[1000];
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
#include <stdlib.h>
#include <math.h>
#include <time.h>

void burn_loop(unsigned int nloops)
{
  long double vf[100][100];
  int vi[1000];
  int i,j,d;

  srand(time(NULL));

  for (i = 0; i < 100; i++)
  {
    d = rand();
    while (d == 0)
      d = rand();

    vf[i][0] = (rand()-d)*rand()/((long double)d);
    vf[0][i] = (d-rand())*rand()/((long double)d);
  }

  for (i = 0; i < 1000; i++)
  {
    vi[i] = rand();
  }


  while (nloops--)
  {
    for (i = 1; i < 100; i++)
      for (j = 1; j < 100; j++)
      {
        int i = 0;
        long double tmp1,
                    tmp2,
                    tmp3,
                    tmp4;

        tmp1 = vi[0];
        do
        {
          tmp2 = (vf[i-1][j] + vf[i][j-1]) / tmp1;
          i = (i+1) % 1000;
          tmp1 += vi[i];
        } while (isinf(tmp2));
        
        do
        {
          tmp3 = y1l(fabsl(vf[i-1][j-1]/tmp1));
          i = (i+1) % 1000;
          tmp1 += vi[i];
        } while (isinf(tmp3));

        tmp4 = cosl(vf[i-1][j])*atanl(vf[i][j-1]) + sinl(vf[i-1][j-1]);

        vf[i][j] = (tmp2 + tmp3)/tmp4;
      }

  }
}

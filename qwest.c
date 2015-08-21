/* rewrite of MULTI5.F */

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>

#define nprimes 97
#define maxn  1000

int primes[] = {  2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
                 31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
                 73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
                127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
                179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
                233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
                283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
                353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
                417, 421, 431, 433, 439, 443, 449, 457, 461, 463,
                467, 479, 487, 491, 499, 503, 509 };

int plist[nprimes];
int otable[nprimes];
int o1list[nprimes];    // list of primes with ord(p,b) = 1
int b = 2;
int nplist = 0;
int opmax = 0;
int o1max = 0;
unsigned int *rpntab;
uint64_t kmin, kmax, kstep;
int low, high;
FILE *zerofile;
FILE *lowfile;
FILE *highfile;

int powmod(int b, int n, int m)    /* powmod = b^n mod m */
{
  int res = 1;
  int d;
  while (n > 0)
  {
    d = n%2;
    if (d == 1)
      res = (b*res)%m;
    b = (b*b)%m;
    n = (n-d)/2;
  }
  return res;
}

int ord(int a, int b)
{
  int k = 1;
  int res = b%a;
  while (res != 1)
  {
    k += 1;
    res *= b;
    res = res%a;
  }
  return k;
}

void init_plist(void)
{
  int i;
  int p;
  int count = 0;
  int ocount = 0;
  int o1count = 0;
  for (i=0; i<nprimes; i++)
  {
    p = primes[i];
    if (p%b != 0)    
    {
      plist[count] = p;
      otable[count] = ord(p,b);
      ocount += otable[count];
//      printf("p = %d otable[%d] = %d\n", p, count, otable[count]);
      if (otable[count] == 1)
        o1list[o1count++] = p;
      count++;
    }
  }
  nplist = count;
  opmax = ocount;
  o1max = o1count;
//  printf("o1max = %d\n", o1max);
}

void init_rpntab(void)
{
  int i, j, k, n, p;
  unsigned int *pntab;
  pntab = (unsigned int *) calloc(opmax+1, sizeof(unsigned int));
  rpntab = (unsigned int *) calloc(opmax+1, sizeof(unsigned int));
  if ((pntab == NULL) || (rpntab == NULL))
  {
    printf("memory allocation error!\n");
    exit(1);
  }

  j = 0;
  for (i=0; i<nplist; i++)
  {
    k = 1;
    p = plist[i];
    for (n=0; n<otable[i]; n++)
    {
//      printf("k = %d p = %d otable(%d) = %d\n", k, p, i, otable[i]);
      k = (b*k)%p;
      j++;
      pntab[j] = k;
//      printf("pntab(%d) = %d\n", j, k);
    }
    rpntab[j] = pntab[j];
    for (n=1; n<otable[i]; n++)
      rpntab[j-n] = pntab[j-otable[i]+n];
  }
  rpntab[0] = 1;
// rpntab[opmax] = 1;
  free(pntab);
}

void open_files(void)
{
  zerofile = fopen("zero.txt", "a");
  lowfile  = fopen("low.txt",  "a");
  highfile = fopen("high.txt", "a");
  if ((!zerofile) || (!lowfile) || (!highfile))
  {
    printf("error creating/opening files!\n");
    exit(1);
  }
}

void close_files(void)
{
  fclose(zerofile);
  fclose(lowfile);
  fclose(highfile);
}

void sieve(void)
{
  uint64_t k;
  int i, j, l, m, n, o, p;
  int count;
  unsigned int ks;
  bool skip = false;
  bool *remain;
  bool *full_remain;
  remain = (bool *) calloc(maxn, sizeof(bool));
  full_remain = (bool *) calloc(maxn, sizeof(bool));
  if ((remain == NULL) || (full_remain == NULL))
  {
    printf("memory allocation error!\n");
    exit(1);
  }
  for (n=0; n<maxn; n++)
    full_remain[n] = true;

  for (k=kmin; k<=kmax; k+=kstep)
  {
    skip = false;
    if ((k%b) == 0)
      skip = true;

    for (i=0; i<o1max; i++)
    {
      p = o1list[i];
      if (((k%p)*(b%p))%p == 1)
      {
        skip = true;
        break;
      }
    }

    if (!skip)
    {
//      for (n=0; n<maxn; n++)
//        remain[n] = true;
      memcpy(remain, full_remain, maxn*sizeof(bool));
      j = 0;
      for (i=0; i<nplist; i++)
      {
        p = plist[i];
        o = otable[i];
        ks = k%p;
        if (ks > 0)
          for (n=1; n<o+1; n++)
          {
            m = n%p;
            if (ks == rpntab[j+m])
            {
              for (l=n; l<=maxn; l+=o)
                remain[l-1] = false;
              break;
            }
          }
        j += o;
      }
      count = 0;
      for (n=0; n<maxn; n++)
        if (remain[n] == true)
          count++;
      if (count == 0)
        fprintf (zerofile, "%20" PRIu64 " %4d\n", k, count);
      else
      {
        if (count <= low)
          fprintf (lowfile,  "%20" PRIu64 " %4d\n", k, count);
        if (count >= high)
          fprintf (highfile, "%20" PRIu64 " %4d\n", k, count);
      }
/*
      for (n=0; n<maxn; n++)
        if (remain[n] == true)
          printf("%d\n", n+1);
*/
    }
  }
  free(remain);
  free(full_remain);
}

int main(int argc, char *argv[])
{
  int option;
  char *ptr;

/* default values */
  b     =       2;
  kmin  =       1;
  kmax  =      10;
  kstep =       2;
  low   =     333;
  high  =     334;

  while ((option = getopt(argc, argv, "b:k:K:s:l:h:")) >= 0)
    switch (option)
    {
      case 'b' : b = atoi(optarg);
                 break;
      case 'k' : kmin = strtoull(optarg, &ptr, 10);
                 break;
      case 'K' : kmax = strtoull(optarg, &ptr, 10);
                 break;
      case 's' : kstep = strtoull(optarg, &ptr, 10);
                 break;
      case 'l' : low = atoi(optarg);
                 break;
      case 'h' : high = atoi(optarg);
                 break;
      case '?' : return 1;
    }

//  printf("b = %d\n", b);
//  printf ("k = %" PRIu64 "-%" PRIu64 "\n", kmin, kmax);
  init_plist();
  init_rpntab();
  open_files();
  sieve();
  close_files();

  return 0;
}

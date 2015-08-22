/* rewrite of MULTI5.F */

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>

int *primes;
int *plist;
int *otable;
int *o1list;    // list of primes with ord(p,b) = 1
int *nmap;      // array of n for which p | k*b^n-1
int b = 2;
int nplist = 0;
int opmax = 0;
int o1max = 0;
int nprimes;
int maxp = 512;
int maxord = 512;
int maxn = 1000;
uint64_t kmin, kmax, kstep;
int low, high;
FILE *zerofile;
FILE *lowfile;
FILE *highfile;

int Erathosthenes(int pmax)
{
  int numbers[pmax];
  int i, j, count;
  for (i=0; i<pmax-1; i++)
    numbers[i] = i+2;
// sieve
  for (i=0; i<pmax-1; i++)
    if (numbers[i] > 0)
      for (j=2*numbers[i]-2; j<pmax; j+=numbers[i])
        numbers[j] = 0;
// count the primes
  count = 0;
  for (i=0; i<pmax-1; i++)
    if (numbers[i] > 0)
      count++;
  primes = (int *) calloc(count, sizeof(int));
// transfer the primes to their own array
  j = 0;
  for (i=0; i<pmax-1; i++)
    if (numbers[i] > 0)
      primes[j++] = numbers[i];
  return j;
}

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
  int o;
  int count = 0;
  int ocount = 0;
  int o1count = 0;
  plist = (int *) calloc(nprimes, sizeof(int));
  otable = (int *) calloc(nprimes, sizeof(int));
  o1list = (int *) calloc(nprimes, sizeof(int));
  for (i=0; i<nprimes; i++)
  {
    p = primes[i];
    if (p%b != 0)    
    {
      o = ord(p,b);
      if ((o > 1) && (o <= maxord))
      {
        plist[count] = p;
        otable[count] = o;
        ocount += o;
        count++;
//      printf("p = %d otable[%d] = %d\n", p, count, otable[count]);
      }
      if (o == 1)
        o1list[o1count++] = p;
    }
  }
  nplist = count;
  opmax = ocount;
  o1max = o1count;
  printf("opmax = %d\n", opmax);
  printf("o1max = %d\n", o1max);
}

void init_nmap(void)
{
  int i, j, k, n, p, m;
  nmap = (int *) calloc(nprimes*maxp, sizeof(int));
  if (nmap == NULL)
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
      m = otable[i]-n-1;
      if (m == 0)
        m+=otable[i];
      nmap[i*maxp+k] = m;
      j++;
    }
  }
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
  int i, j, l, m, n, o, p, om, nm;
  int count;
  int kmodp[nplist];    // precomputed k%p for all p in plist
  int kbmodp[o1max];    // precomputed (k*b^1)%p for all p in o1list
  int bksmodp[o1max];   // precomputed (b*kstep)%p for all p in o1list
  int kmodb;            // precomputed k%b
  int ks;
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

  kmodb = (kmin-kstep)%b;
  for (i=0; i<nplist; i++)
  {
    p = plist[i];
    kmodp[i] = (kmin-kstep)%p;
  }
  for (i=0; i<o1max; i++)
  {
    p = o1list[i];
    kbmodp[i] = (b*(kmin-kstep)%p)%p;
    bksmodp[i] = (b*kstep)%p;
  }
  
  for (k=kmin; k<=kmax; k+=kstep)
  {
    for (i=0; i<nplist; i++)
    {
      kmodp[i] += kstep;
      while (kmodp[i] >= plist[i])
        kmodp[i] -= plist[i];
    }

    for (i=0; i<o1max; i++)
    {
      kbmodp[i] += bksmodp[i];
      while (kbmodp[i] >= o1list[i])
        kbmodp[i] -= o1list[i];
    }

    skip = false;
    kmodb += kstep;
    while (kmodb >= b)
      kmodb -= b;
    if (kmodb == 0)
      skip = true;

    if (!skip)
    {
      for (i=0; i<o1max; i++)
      {
        if (kbmodp[i] == 1)
        {
          skip = true;
          break;
        }
      }
    }

    if (!skip)
    {
//      for (n=0; n<maxn; n++)
//        remain[n] = true;
      memcpy(remain, full_remain, maxn*sizeof(bool));
      for (i=0; i<nplist; i++)
      {
        p = plist[i];
        o = otable[i];
        ks = kmodp[i];
        if (ks > 0)
        {
          n = nmap[i*maxp+ks];
          if (n > 0)
            for (l=n; l<=maxn; l+=o)
              remain[l-1] = false;
        }
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
  b      =       2;
  kmin   =       1;
  kmax   =      10;
  kstep  =       2;
  low    =     333;
  high   =     334;
  maxp   =     512;
  maxord =     512;
  maxn   =    1000;

  while ((option = getopt(argc, argv, "b:k:K:s:l:h:p:o:n:")) >= 0)
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
      case 'p' : maxp = atoi(optarg);
                 break;
      case 'o' : maxord = atoi(optarg);
                 break;
      case 'n' : maxn = atoi(optarg);
                 break;
      case '?' : return 1;
    }

//  printf("b = %d\n", b);
//  printf ("k = %" PRIu64 "-%" PRIu64 "\n", kmin, kmax);
  nprimes = Erathosthenes(maxp);
  printf("no. of primes = %d\n", nprimes);
  printf("largest prime = %d\n", primes[nprimes-1]);
  maxp = primes[nprimes-1];

/*
  int i;
  for (i=0; i<nprimes; i++)
    printf("%4d", primes[i]);
  printf("\n");
*/
   
  init_plist();
  init_nmap();

/*
  int i, j;
  for (i=0; i<nprimes; i++)
  {
    for (j=0; j<maxp; j++)
      printf("%3d", nmap[i*maxp+j]);
    printf("\n");
  }
*/

  open_files();
  sieve();
  close_files();

  return 0;
}

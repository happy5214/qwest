/* rewrite of MULTI5.F */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <signal.h>

#include "int128.h"

#define REPORT_INTERVAL 5000000

#define BUFFER_SIZE 50

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
int slice = 0;
int modulus = 1;
bool quiet = false;
bool stop = false;
bool slicing = false;
bool riesel = false;
bool ignore_zeros = false;
uint128_t kmin, kmax, kstep;
int low, high;
FILE *zerofile;
FILE *lowfile;
FILE *highfile;

char buffer[BUFFER_SIZE];

void terminate(int signum)
{
  signal (signum, SIG_IGN);
  stop = true;
}

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

void init_plist(bool skip_kstep_factors)
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
    if ((p%b != 0) && (b%p != 0) && ((skip_kstep_factors && (kstep%p != 0)) || !skip_kstep_factors))
    {
      o = ord(p,b);
      if ((o > 1) && (o <= maxord))
      {
        plist[count] = p;
        otable[count] = o;
        ocount += o;
//        printf("p = %d otable[%d] = %d\n", p, count, otable[count]);
        count++;
      }
      if (o == 1)
        o1list[o1count++] = p;
    }
  }
  nplist = count;
  opmax = ocount;
  o1max = o1count;
  if (!quiet)
  {
    printf("opmax = %d\n", opmax);
    printf("o1max = %d\n", o1max);
  }
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
  if (!ignore_zeros)
    zerofile = fopen("zero.txt", "a");
  lowfile  = fopen("low.txt",  "a");
  highfile = fopen("high.txt", "a");
  if ((!ignore_zeros && !zerofile) || (!lowfile) || (!highfile))
  {
    printf("error creating/opening files!\n");
    exit(1);
  }
}

void read_checkpoint()
{
  FILE *file;
  uint128_t k;
  if ((file = fopen("checkpoint.txt", "r")) == NULL)
    return;

  if (fscanf(file,"%s", buffer) == 1)
  {
    k = strtou128(buffer, NULL, 10);
    if (kmin < k)
    {
      kmin = k;
      snprint_u128(buffer, BUFFER_SIZE, k);
      printf("Resuming from checkpoint k = %s\n", buffer);
    }
  }
  fclose(file);
  remove("checkpoint.txt");
}

void write_checkpoint(uint128_t k)
{
  FILE *file;
  file = fopen("checkpoint.txt", "w");
  if (!file)
  {
    printf("error creating checkpoint file!\n");
    exit(1);
  }
  snprint_u128(buffer, BUFFER_SIZE, k);
  fprintf(file, "%s\n", buffer);
  fclose(file);
}

void close_files(void)
{
  if (!ignore_zeros)
    fclose(zerofile);
  fclose(lowfile);
  fclose(highfile);
}

void sieve(void)
{
//  uint64_t k;
  uint128_t k;

  double to_percent;    // for percentage calculation
  int i, l, n, o, p;
  int countdown;        // for progress indicator
  int count;
  int kmodp[nplist];    // precomputed k%p for all p in plist
  int kstepmodp[nplist];   // precomputed kstep%p for all p in plist
  int kbmodp[o1max];    // precomputed (k*b^1)%p for all p in o1list
  int bksmodp[o1max];   // precomputed (b*kstep)%p for all p in o1list
  int kmodb;            // precomputed k%b
  int kstepmodb;        // precomputed kstep%b
  int ks;
  int pos;
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

//  kmodb = ((kmin-kstep)%b+b)%b;
  kmodb = (b + kmin%b - kstep%b)%b;
  kstepmodb = kstep%b;
  if (!riesel)
  {
    kmodb = b - kmodb;
    kstepmodb = b - kstepmodb;
  }

  for (i=0; i<nplist; i++)
  {
    p = plist[i];
//    kmodp[i] = (p+(kmin-kstep)%p)%p;
//    kmodp[i] = (p + kmin%p - kstep%p)%p;
    if (riesel)
      kmodp[i] = (p + kmin%p - kstep%p)%p;
    else
      kmodp[i] = p - (p + kmin%p - kstep%p)%p;
//    printf("p = %d, kmodp = %d\n", p, kmodp[i]);
//    kstepmodp[i] = kstep%p;
    if (riesel)
      kstepmodp[i] = kstep%p;
    else
      kstepmodp[i] = p - kstep%p;
  }
  for (i=0; i<o1max; i++)
  {
    p = o1list[i];
//    kbmodp[i] = (b*((kmin-kstep)%p+p))%p;
//    kbmodp[i] = (b*(p + kmin%p - kstep%p))%p;
//    bksmodp[i] = (b*kstep)%p;
    if (riesel)
    {
      kbmodp[i] = (b*(p + kmin%p - kstep%p))%p;
      bksmodp[i] = (b*kstep)%p;
    }
    else
    {
      kbmodp[i] = p - (b*(p + kmin%p - kstep%p))%p;
      bksmodp[i] = p - (b*kstep)%p;
    }
//    printf("p = %d, kbmodp = %d, bksmodp = %d\n", p, kbmodp[i], bksmodp[i]);
  }
  
  to_percent = 100.0/(double)(kmax - kmin);
// adjust kmax accordingly so that: kmax = kmin + x*kstep
  kmax -= (kmax-kmin)%kstep;
  if (kmax + kstep < kmax)    // to prevent overflow at 2^128-1
  {
    kmax -= kstep;
    snprint_u128(buffer, BUFFER_SIZE, kmax);
    printf ("kmax (adjusted) = %s\n", buffer);
  }

  countdown = REPORT_INTERVAL;

  for (k=kmin; k<=kmax; k+=kstep)
  {
    for (i=0; i<nplist; i++)
    {
      kmodp[i] += kstepmodp[i];
      if (kmodp[i] >= plist[i])
        kmodp[i] -= plist[i];
    }

    for (i=0; i<o1max; i++)
    {
      kbmodp[i] += bksmodp[i];
      if (kbmodp[i] >= o1list[i])
        kbmodp[i] -= o1list[i];
    }

    skip = false;
    kmodb += kstepmodb;
    if (kmodb >= b)
      kmodb -= b;
    if (kmodb == 0)
      skip = true;

    if (!skip)
    {
      for (i=0; i<o1max; i++)
      {
        if (kbmodp[i] == 1)     // Riesel
//        if (kbmodp[i] == plist[i]-1)  // original Sierpinski (no inverting)
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
      pos = 0;
      for (i=0; i<nplist; i++)
      {
        p = plist[i];
        o = otable[i];
        ks = kmodp[i];
        if (ks > 0)
        {
          n = nmap[pos+ks];
//          n = nmap[i*maxp+ks];
          if (n > 0)
            for (l=n; l<=maxn; l+=o)
              remain[l-1] = false;
        }
        pos += maxp;
      }
      count = 0;
      if (slicing)
      {
        for (n=slice; n<maxn; n+=modulus)
          if (remain[n] == true)
            count++;
      }
      else
      {
        for (n=0; n<maxn; n++)
          if (remain[n] == true)
            count++;
      }

      if (count == 0)
      { 
        if (!ignore_zeros)
        {
          n = snprint_u128(buffer, BUFFER_SIZE, k);
          fprintf (zerofile, "%40s %4d\n", buffer, count);
        }
      }
      else
      {
        if (count <= low)
        {
          n = snprint_u128(buffer, BUFFER_SIZE, k);
          fprintf (lowfile,  "%40s %4d\n", buffer, count);
        }
        if (count >= high)
        {
          n = snprint_u128(buffer, BUFFER_SIZE, k);
          fprintf (highfile, "%40s %4d\n", buffer, count);
        }
      }

/*
      for (n=0; n<maxn; n++)
        if (remain[n] == true)
          printf("%d\n", n+1);
*/
    }
    if (stop)
    {
      n = snprint_u128(buffer, BUFFER_SIZE, k);
      printf("Terminating at k = %s\n", buffer);
      write_checkpoint(k+kstep);
      break;
    }
    countdown--;
    if ((countdown == 0) && !quiet)
    {
      n = snprint_u128(buffer, BUFFER_SIZE, k);
      printf("Tested up to k = %s (%.2f%% done)\n", buffer, (k-kmin)*to_percent);
      countdown = REPORT_INTERVAL;
    }
  }
  free(remain);
  free(full_remain);
}

int main(int argc, char *argv[])
{
  int option;
//  char *ptr;

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
  quiet  =   false;
  stop   =   false;
  slicing =  false;
  riesel =   false;
  slice  =       0;
  modulus =      1;

  while ((option = getopt(argc, argv, "b:k:K:s:l:h:p:o:n:e:m:qzr")) >= 0)
    switch (option)
    {
      case 'b' : b = atoi(optarg);
                 break;
      case 'k' : kmin = strtou128(optarg, NULL, 10);
// kmin = strtoull(optarg, &ptr, 10);
                 break;
      case 'K' : kmax = strtou128(optarg, NULL, 10);
// kmax = strtoull(optarg, &ptr, 10);
                 break;
      case 's' : kstep = strtou128(optarg, NULL, 10);
// kstep = strtoull(optarg, &ptr, 10);
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
      case 'e' : slice = atoi(optarg)-1;
                 slicing = true;
                 break;
      case 'm' : modulus = atoi(optarg);
                 slicing = true;
                 break;
      case 'q' : quiet = true;
                 break;
      case 'z' : ignore_zeros = true;
                 break;
      case 'r' : riesel = true;
                 break;
      case '?' : return 1;
    }

  read_checkpoint();
  signal(SIGINT, terminate);
  signal(SIGTERM, terminate);
//  printf("b = %d\n", b);
//  snprint_u128(buffer, BUFFER_SIZE, kmin);
//  printf ("k = %s", buffer);
//  snprint_u128(buffer, BUFFER_SIZE, kmax);
//  printf ("-%s", buffer);
  nprimes = Erathosthenes(maxp);
  maxp = primes[nprimes-1];
  if (!quiet)
  {
    printf("no. of primes = %d\n", nprimes);
    printf("largest prime = %d\n", maxp);
  }

/*
  int i;
  for (i=0; i<nprimes; i++)
    printf("%4d", primes[i]);
  printf("\n");
*/
   
  init_plist(slicing);
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

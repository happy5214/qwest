# qwest - Quick weight estimator

This program is a quick estimator of prime density for Riesel and Sierpiński *k*
values. It is intended to search large ranges at once for *k*'s with weights
above or below a certain limit, with the goal of further weight testing with
the [Nash weight calculator](https://github.com/happy5214/nash) and eventual
primality testing. Unlike the Nash weight calculator, this program uses an
ordinary sieve to achieve a much faster runtime.

qwest is based on the multi5 Fortran program originally released on the
[MersenneForum](https://www.mersenneforum.org/showpost.php?p=99532) in 2007
by Thomas Ritschel. He is responsible for most of this rewrite's code as well.
This program is released under the MIT license.

# Compiling on a Linux system

To build on a Linux system, simply type `make`.

# Usage example

```
./qwest -r -b 2 -k 1 -K 1000000 -s 2 -n 1000 -p 12000 -o 256 -l 15 -h 650 -q
```

# Explanation of the command line arguments

- `-r`: Riesel side (*k*\*2^n-1), if not given: Sierpiński side (*k*\*2^n+1)
- `-b`: *base* (default: 2)
- `-k`: *kmin* (default: 1)
- `-K`: *kmax* (default: 10)
- `-s`: *kstep* (step size for *k*)      (default: 2)
- `-n`: *nmax* (*n* range is 1...*nmax*)   (default: 1000)
- `-p`: *pmax* (use all primes up to *pmax* for the sieve)
- `-o`: *omax* (use only primes with multiplicative order ord(*b*,*p*) ≤ *omax*)
- `-l`: *wlow* (print only candidates with weight ≤ *wlow*)
- `-h`: *whigh* (print only candidates with weight ≥ *whigh*)
- `-q`: quiet (do not print start-up debug information)
- `-z`: skip printing the zero weight candidates

There will be no screen output (only some information about the primes used
for the sieve on start-up if the quiet (`-q`) option isn't set.

The results are written into the following 3 files:

- `low.txt`   (containing the low weight candidates)
- `high.txt`  (containing the high weight candidates)
- `zero.txt`  (containing the zero weight candidates)

The amount of output can be controlled by the `-l` and `-h` options.
To skip storing the low weights entirely one can set *wlow* to 0: `-l 0`
To skip the high weights one can set *whigh* to *nmax* (or higher): `-h 1000`

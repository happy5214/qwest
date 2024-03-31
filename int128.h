#ifndef INT128_H
#define INT128_H 1

#include <stdio.h>

typedef unsigned __int128 uint128_t;

uint128_t strtou128(const char *s);
int snprint_u128(char *buffer, size_t size, uint128_t u128);

#endif // INT128_H

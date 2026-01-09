#ifndef INT128_H
#define INT128_H 1

#include <iostream>

typedef __uint128_t uint128_t;

#define UINT128_MAX ((uint128_t) 0) - 1

uint128_t strtou128(const char * const string, const char **end, const int base);
int snprint_u128(char * const buffer, const size_t size, const uint128_t u128);

std::ostream& operator<<(std::ostream& out, __uint128_t val);
std::ostream& operator<<(std::ostream& out, __int128_t val);

#endif // INT128_H

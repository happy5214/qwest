#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

#include "int128.h"

#define P10_UINT64 10000000000000000000ULL

uint128_t strtou128(const char *string) {
  const int base = 10;
  uint128_t number = 0;

  for(; *string; ++string) {
    unsigned char digit = *string;
    digit -= '0';
    number *= (uint128_t) base;
    number += (uint128_t) digit;
//    sprint_u128(buffer, number); 
//    printf("%d %s\n", c, buffer);
  }

  return number;
}

int snprint_u128(char *buffer, size_t size, uint128_t u128) {
  int charactersWritten = 0;

  if (u128 > UINT64_MAX) {
    const uint128_t leading = u128 / P10_UINT64;
    const uint64_t trailing = u128 % P10_UINT64;

    if (leading > UINT64_MAX) {
      const uint64_t u64Leading = leading / P10_UINT64;
      const uint64_t u64Trailing = leading % P10_UINT64;
      charactersWritten = snprintf(buffer, size, "%" PRIu64 "%.19" PRIu64 "%.19" PRIu64, u64Leading, u64Trailing, trailing);
    } else {
      const uint64_t u64 = leading;
      charactersWritten = snprintf(buffer, size, "%" PRIu64 "%.19" PRIu64, u64, trailing);
    }
  } else {
    const uint64_t u64 = u128;
    charactersWritten = snprintf(buffer, size, "%" PRIu64, u64);
  }

  return charactersWritten;
}

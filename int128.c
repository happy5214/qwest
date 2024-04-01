#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

#include "int128.h"

#define P10_UINT64 10000000000000000000ULL

uint128_t strtou128(const char * const string, const char **end, int base)
{
  uint128_t number = 0;
  const char *stringPointer = string;

  if (!base) {
    base = 10;
  }

#ifdef DEBUG
  char buffer[50];
  int c = 0;
#endif

  for(; *stringPointer; ++stringPointer)
  {
    unsigned char digit = *stringPointer;
    if (isspace(digit))
      continue;

    digit -= '0';

    if (digit >= 0 && digit <= 9)
    {
      const uint128_t nextNumber = number * (uint128_t) base + (uint128_t) digit;
      if (nextNumber < number)
      {
        number = UINT128_MAX;
        errno = ERANGE;
        break;
      }
      number = nextNumber;
    }
    else
    {
      if (end)
        *end = stringPointer;
      break;
    }

#ifdef DEBUG
    snprint_u128(buffer, 50, number);
    printf("%d %s\n", c, buffer);
    c++;
#endif
  }

  return number;
}

int snprint_u128(char * const buffer, const size_t size, const uint128_t u128)
{
  int charactersWritten = 0;

  if (u128 > UINT64_MAX)
  {
    const uint128_t leading = u128 / P10_UINT64;
    const uint64_t trailing = u128 % P10_UINT64;

    if (leading > UINT64_MAX)
    {
      const uint64_t u64Leading = leading / P10_UINT64;
      const uint64_t u64Trailing = leading % P10_UINT64;
      charactersWritten = snprintf(buffer, size, "%" PRIu64 "%.19" PRIu64 "%.19" PRIu64, u64Leading, u64Trailing, trailing);
    }
    else
    {
      const uint64_t u64 = leading;
      charactersWritten = snprintf(buffer, size, "%" PRIu64 "%.19" PRIu64, u64, trailing);
    }
  }
  else
  {
    const uint64_t u64 = u128;
    charactersWritten = snprintf(buffer, size, "%" PRIu64, u64);
  }

  return charactersWritten;
}

#include <cctype>
#include <cerrno>
#include <cinttypes>
#include <cstdint>
#include <cstdio>

#include "int128.h"

#define P10_UINT64 10000000000000000000ULL

uint128_t strtou128(const char * const string, const char **end, const int base)
{
#ifdef DEBUG
  char buffer[50];
  int charIndex = 0;
#endif

  uint128_t number = 0;
  const char *stringPointer = string;

  if (base != 10)
  {
    errno = EINVAL;
    return 0;
  }

  while (isspace(*stringPointer))
    stringPointer++;

  for(; *stringPointer; ++stringPointer)
  {
    unsigned char digit = *stringPointer;

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
    printf("%d %s\n", charIndex, buffer);
    charIndex++;
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

// Source - https://stackoverflow.com/a/72651639
// Posted by Mark Adler, modified by community. See post 'Timeline' for change history
// Retrieved 2026-01-09, License - CC BY-SA 4.0

#include <iostream>
#include <iomanip>

// Write the 128-bit integer val to out, with a minus sign if decimal and neg
// is true. Obey all of the ostream settings of out for integer display: octal
// or hexadecimal, upper case letters, plus sign, fill character and width, and
// fill placement.
static void out128(std::ostream& out, __uint128_t val, int neg) {
    // Note if the number is zero. (No hex or octal prefix in this case.)
    auto zero = val == 0;

    // Note if upper-case letters requested.
    auto state = out.flags();
    auto upper = (state & std::ios_base::uppercase) != 0;

    // Set base for digits.
    unsigned base = state & std::ios_base::hex ? 16 :
                    state & std::ios_base::oct ? 8 :
                    10;

    // Space for digits and prefix. Generate digits starting at the end of the
    // string, going backwards. num will be the digit string. Terminate it.
    char str[47];
    auto end = str + sizeof(str), num = end;
    *--num = 0;

    // Compute and place digits in base base.
    do {
        char dig = val % base;
        val /= base;
        dig += dig < 10 ? '0' : (upper ? 'A' : 'a') - 10;
        *--num = dig;
    } while (val);

    // Prepend octal number with a zero if requested.
    if (state & std::ios_base::showbase && base == 8 && !zero)
        *--num = '0';

    // pre will be the prefix string. Terminate it.
    auto pre = num;
    *--pre = 0;

    // Put a plus or minus sign in the prefix as appropriate.
    if (base == 10) {
        if (neg)
            *--pre = '-';
        else if (state & std::ios_base::showpos)
            *--pre = '+';
    }

    // Prefix a hexadecimal number if requested.
    else if (state & std::ios_base::showbase && base == 16 && !zero) {
        *--pre = upper ? 'X' : 'x';
        *--pre = '0';
    }

    // Compute the number of pad characters and get the fill character.
    auto len = (num - pre) + (end - num) - 2;
    auto pad = out.width();
    out.width(0);
    pad = pad > len ? pad - len : 0;
    char fill = out.fill();

    // Put the padding before prefix if neither left nor internal requested.
    if (!(state & (std::ios_base::internal | std::ios_base::left)))
        while (pad) {
            out << fill;
            pad--;
        }

    // Write prefix.
    out << pre;

    // Put the padding between the prefix and the digits if requested.
    if (state & std::ios_base::internal)
        while (pad) {
            out << fill;
            pad--;
        }

    // Write digits.
    out << num;

    // Put number to the left of padding, if requested.
    if (state & std::ios_base::left)
        while (pad) {
            out << fill;
            pad--;
        }
}

// Overload << for an unsigned 128-bit integer.
std::ostream& operator<<(std::ostream& out, __uint128_t val) {
    out128(out, val, 0);
    return out;
}

// Overload << for a signed 128-bit integer. Negation of the most negative
// signed value gives the correct unsigned absolute value.
std::ostream& operator<<(std::ostream& out, __int128_t val) {
    auto state = out.flags();
    if (val < 0 && !(state & (std::ios_base::hex | std::ios_base::oct)))
        out128(out, -(__uint128_t)val, 1);
    else
        out128(out, (__uint128_t)val, 0);
    return out;
}

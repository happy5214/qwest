#include <vector>

#include "qwest.h"

std::vector<int> Erathosthenes(const int p_max) {
	std::vector<int> primes;
	std::vector<bool> numbers(p_max, true);

	// Sieve
	for (int i = 0; i < p_max - 1; i++) {
		if (numbers[i]) {
			for (int j = 2 * i + 2; j < p_max; j += (i + 2)) {
				numbers[j] = false;
			}
			primes.push_back(i + 2);
		}
	}
	return primes;
}

/* powmod = b^n mod m */
int powmod(int base, int exponent, const int modulus) {
	int result = 1;
	int low_bit;
	while (exponent > 0) {
		low_bit = exponent & 2;
		if (low_bit == 1) {
			result = (base * result) % modulus;
		}
		base = (base * base) % modulus;
		exponent >>= 1;
	}
	return result;
}

int ord(const int a, const int b) {
	int k = 1;
	int result = b % a;
	while (result != 1) {
		k += 1;
		result *= b;
		result = result % a;
	}
	return k;
}

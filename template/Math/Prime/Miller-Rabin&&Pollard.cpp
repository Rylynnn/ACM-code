#include<cstdio>
#include<typeinfo>
#include<cstdlib>
#include<algorithm>
using namespace std;
typedef unsigned long long ULL;
typedef unsigned int UInt;
const UInt base1[] = { 2, 7, 61, 0 };
const UInt base2[] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022, 0 };
const UInt prime[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53 };
template <typename T>
inline T add(T a, T b, T mod){
	return a + b - (a + b >= mod ? mod : 0);
}
inline UInt mul(UInt a, UInt b, UInt mod){
	return (ULL)a * b % mod;
}
ULL mul(ULL a, ULL b, ULL mod){
	ULL ret = 0;
	for (ULL t = a; b; b >>= 1){
		if (b & 1)ret = add(ret, t, mod);
		t <<= 1;
		if (t >= mod)t -= mod;
	}
	return ret;
}
template <typename T>
T power(T a, T b, T mod){
	T ret = 1;
	for (T t = a; b; b >>= 1){
		if (b & 1)ret = mul(ret, t, mod);
		t = mul(t, t, mod);
	}
	return ret;
}
//n为小于2^63的非1奇数，正确性100%
template <typename T>
bool millerRabin(T n)
{
	int s = 0;
	T r = n;
	for (r--; !(r & 1); r >>= 1)s++;
	for (const UInt *base = typeid(T) == typeid(UInt) ? base1 : base2; *base; base++){
		T t = power(*base % n, r, n);
		if (t == 0 || t == 1 || t == n - 1)continue;
		for (int j = 1; j < s; j++){
			t = mul(t, t, n);
			if (t == 1)return false;
			if (t == n - 1)break;
		}
		if (t != n - 1)return false;
	}
	return true;
}
template <typename T>
bool checkPrime(T n)
{
	if (n == 1)return false;
	for (int i = 0; i < sizeof(prime) / sizeof(int); i++){
		if (n % prime[i] == 0)return n == prime[i];
	}
	return millerRabin(n);
}
template <typename T>
T gcd(T x, T y){
	return y ? gcd(y, x % y) : x;
}
template <typename T>
T pollard(T n)
{
	if (millerRabin(n))return n;
	while (1){
		T x = rand() % n, y = x, c = rand() % (n - 1) + 1;
		for (UInt i = 1, j = 2;; i++){
			if (i == j){ j *= 2; y = x; }
			x = add(mul(x, x, n), c, n);
			T d = gcd(x - y + n, n);
			if (d != 1){
				if (d != n)return d;
				break;
			}
		}
	}
}
ULL factor[64];
int factorNum;
void calFactorInternal(ULL n)
{
	ULL d;
	d = n >> 32 ? pollard(n) : pollard((UInt)n);
	if (d == n){ factor[factorNum++] = d; return; }
	calFactorInternal(d);
	calFactorInternal(n / d);
}
void calFactor(ULL n)
{
	factorNum = 0;
	for (int i = 0; i < sizeof(prime) / sizeof(int); i++){
		while (n % prime[i] == 0){
			n /= prime[i];
			factor[factorNum++] = prime[i];
		}
	}
	if (n != 1)calFactorInternal(n);
	sort(factor, factor + factorNum);
}
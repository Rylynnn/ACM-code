#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
inline int mul(int a, int b, int mod){ return (long long)a*b%mod; }
int power(int a, int b, int mod){
	int ret = 1;
	for (int t = a; b; b >>= 1){
		if (b & 1)ret = mul(ret, t, mod);
		t = mul(t, t, mod);
	}
	return ret;
}
int cal_root(int mod)
{
	int factor[20], num = 0, m = mod - 1, s = m;
	for (int i = 2; i * i <= s; i++){
		if (s % i == 0){
			factor[num++] = i;
			while (s % i == 0)s /= i;
		}
	}
	if (s != 1)factor[num++] = s;
	for (int i = 2;; i++){
		int j = 0;
		for (; j < num && power(i, m / factor[j], mod) != 1; j++);
		if (j == num)return i;
	}
}
template<int MOD, int ROOT>
void fft_main(int a[], int len, bool reverse)
{
	for (int i = 1, j = len / 2; i < len - 1; i++) {
		if (i < j) swap(a[i], a[j]);
		for (int k = len; j < k; k >>= 1, j ^= k);
	}
	for (int i = 1, s = 1; s < len; i++, s <<= 1){
		int t = (MOD - 1) / (s * 2);
		int step = power(ROOT, reverse ? MOD - 1 - t : t, MOD);
		for (int j = 0; j < len; j += 2 * s){
			int cur = 1;
			for (int k = j; k < j + s; k++){
				int u = a[k], t = mul(cur, a[k + s], MOD);
				a[k] = (unsigned int)(u + t) % MOD;
				a[k + s] = (unsigned int)(u - t + MOD) % MOD;
				cur = mul(cur, step, MOD);
			}
		}
	}
	if (reverse){
		int t = power(len, MOD - 2, MOD);
		for (int i = 0; i < len; i++)
			a[i] = mul(a[i], t, MOD);
	}
}
//确保数组中的数小于mod(mod<2^30)，数组需留足2^(logn向上取整+1)的空间
//并且mod为形如m*2^k+1的素数，2^k>=2*n
template<int MOD, int ROOT>
void fft(int a[], int b[], int n)
{
	int len = 1;
	while (len < 2 * n)len <<= 1;
	memset(a + n, 0, sizeof(int)*(len - n));
	memset(b + n, 0, sizeof(int)*(len - n));
	fft_main<MOD, ROOT>(a, len, 0);
	fft_main<MOD, ROOT>(b, len, 0);
	for (int i = 0; i < len; i++)
		a[i] = mul(a[i], b[i], MOD);
	fft_main<MOD, ROOT>(a, len, 1);
}
template<int MOD, int ROOT>
void fft(int a[], int n)
{
	int len = 1;
	while (len < 2 * n)len <<= 1;
	memset(a + n, 0, sizeof(int)*(len - n));
	fft_main<MOD, ROOT>(a, len, 0);
	for (int i = 0; i < len; i++)
		a[i] = mul(a[i], a[i], MOD);
	fft_main<MOD, ROOT>(a, len, 1);
}
//附满足条件大整数：167772161, 469762049, 754974721, 2146959361
//原根：3 3 11 19
#include<cstdio>
#include<cmath>
#include<algorithm>
#include<vector>
#include<cstring>
using namespace std;
int len, bit;
int MOD, w[2][32];
inline int add(int a, int b){
	return a + b - (a + b >= MOD ? MOD : 0);
}
inline int sub(int a, int b){
	return a - b + (a - b < 0 ? MOD : 0);
}
inline int mul(int a, int b){
	return (long long)a * b % MOD;
}
int power(int a, int b){
	int ret = 1;
	for (int t = a; b; b >>= 1){
		if (b & 1)ret = mul(ret, t);
		t = mul(t, t);
	}
	return ret;
}
int cal_root(int mod)
{
	for (int i = 2;; i++){
		if (power(i, (mod - 1) / 2) == mod - 1)
			return i;
	}
}
void fft_init(int n, int mod)
{
	MOD = mod;
	bit = (int)log2(n - 0.5) + 2;
	len = 1 << bit;
	w[0][0] = power(cal_root(mod), (mod - 1) / len);
	int i;
	for (i = 1; i < bit; i++)
		w[0][i] = mul(w[0][i - 1], w[0][i - 1]);
	i--;
	w[1][i] = w[0][i];
	for (i--; i >= 0; i--)
		w[1][i] = mul(w[1][i + 1], w[0][i]);
}
void bitReverse(int a[]) {
	for (int i = 1, j = len / 2; i < len - 1; i++) {
		if (i < j) swap(a[i], a[j]);
		int k = len / 2;
		while (j >= k) { j -= k; k >>= 1; }
		if (j < k) j += k;
	}
}
void fft_main(int a[], bool reverse)
{
	bitReverse(a);
	for (int i = 1, s = 1; s < len; i++, s <<= 1){
		int step = w[reverse][bit - i];
		for (int j = 0; j < len; j += 2 * s){
			int cur = 1;
			for (int k = j; k < j + s; k++){
				int u = a[k], t = mul(cur, a[k + s]);
				a[k] = add(u, t);
				a[k + s] = sub(u, t);
				cur = mul(cur, step);
			}
		}
	}
	if (reverse){
		int t = power(len, MOD - 2);
		for (int i = 0; i < len; i++)
			a[i] = mul(a[i], t);
	}
}
//确保数组中的数小于mod(mod<2^30)，数组需留足2^(logn向上取整+1)的空间，后面填充0
//并且mod为形如m*2^k+1的素数，2^k>=2*n
void fft(int a[], int b[], int n, int mod)
{
	fft_init(n, mod);
	fft_main(a, 0); fft_main(b, 0);
	for (int i = 0; i < len; i++)
		a[i] = mul(a[i], b[i]);
	fft_main(a, 1);
}
//确保mod两两互质，retmod任意
void chineseRemainder(const int mod[], int *a[], int ret[], int num, int n, int retMod)
{
	int kk[30], mulMod[30][30], mulModr[30], mulretMod[30];
	for (int i = 0; i < num; i++){
		MOD = mod[i]; mulMod[i][0] = 1;
		for (int j = 1; j <= i; j++)
			mulMod[i][j] = mul(mulMod[i][j - 1], mod[j - 1]);
		mulModr[i] = power(mulMod[i][i], MOD - 2);
	}
	mulretMod[0] = 1; MOD = retMod;
	for (int i = 1; i < num; i++)
		mulretMod[i] = mul(mulretMod[i - 1], mod[i - 1]);
	for (int i = 0; i < n; i++){
		for (int j = 1; j < num; j++){
			MOD = mod[j];
			int sum = a[0][i] % MOD;
			for (int k = 1; k < j; k++)
				sum = add(sum, mul(mulMod[j][k], kk[k]));
			kk[j] = mul(sub(a[j][i] % MOD, sum), mulModr[j]);
		}
		MOD = retMod;
		ret[i] = a[0][i] % MOD;
		for (int j = 1; j < num; j++)
			ret[i] = add(ret[i], mul(kk[j] % MOD, mulretMod[j]));
	}
}
//附满足条件大整数：167772161, 469762049, 754974721

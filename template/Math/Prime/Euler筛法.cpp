#include<cstdio>
#define MAXN 10000001
int minFactor[MAXN];
int prime[2000000], primeNum;
int phi[MAXN];
void calPrime()
{
	for (int i = 2; i < MAXN; i++){
		if (!minFactor[i]){
			prime[primeNum++] = i;
			minFactor[i] = primeNum;
		}
		for (int j = 1; j <= minFactor[i]; j++){
			int t = i * prime[j - 1];
			if (t >= MAXN)break;
			minFactor[t] = j;
		}
	}
}
void calPhi()
{
	phi[1] = 1;
	for (int i = 2; i < MAXN; i++){
		if (!minFactor[i]){
			prime[primeNum++] = i;
			minFactor[i] = primeNum;
			phi[i] = i - 1;
		}
		for (int j = 1;; j++){
			int t = i * prime[j - 1];
			if (t >= MAXN)break;
			minFactor[t] = j;
			if (j == minFactor[i]){
				phi[t] = phi[i] * prime[j - 1];
				break;
			}
			phi[t] = phi[i] * (prime[j - 1] - 1);
		}
	}
}
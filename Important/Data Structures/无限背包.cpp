#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
int n;
int w[1001], c[1001], tmp[30];
long long f[30][2001], g[3003], s;
int main()
{
	scanf("%d%I64d", &n, &s);
	for (int i = 0; i < n; i++)
		scanf("%d%d", &w[i], &c[i]);
	int maxW = *max_element(w, w + n), best = 0, cnt = 0;
	for (int i = 1; i < n; i++){
		if ((long long)c[best] * w[i]<(long long)c[i] * w[best])
			best = i;
	}
	long long num = max(s / w[best] - maxW, 0LL);
	for (int j = s - num * w[best]; j > 0; j = (j - maxW) / 2, cnt++)
		tmp[cnt + 1] = j;
	for (int i = 0; i < n; i++){
		for (int j = w[i]; j <= 3 * maxW + 1; j++)
			g[j] = max(g[j], g[j - w[i]] + c[i]);
	}
	memcpy(f[cnt], g + tmp[cnt], (maxW * 2 + 1)*sizeof(long long));
	while (--cnt > 0){
		int t = tmp[cnt] - 2 * tmp[cnt + 1];
		for (int i = 0; i <= 2 * maxW; i++){
			for (int j = (i + 1) / 2; j <= (i + 1) / 2 + maxW / 2 + 1; j++)
				f[cnt][i] = max(f[cnt][i], f[cnt + 1][j] + f[cnt + 1][t + i - j]);
		}
	}
	printf("%I64d", f[1][0] + num * c[best]);
}
#include<algorithm>
#include<cstdio>
#include<cstring>
using namespace std;
#define MAXN 2000001
#define CHAR 256
int bucket[MAXN];
char s[MAXN];
//#define DC3
#ifndef DC3
int r[2][MAXN];
int tmp[MAXN], sa[MAXN];//存储第i小的字符位置
int *rk;//存储位置i上的字符是第几小的
inline bool cmp(int t1, int t2){
	return s[t1] < s[t2];
}
int suffixArray()
{
	int len = 0, m;
	int *r1 = r[0], *r2 = r[1];
	do sa[len] = len;
	while (s[len++]);
	sort(sa, sa + len, cmp);
	r[0][len - 1] = m = 0;
	for (int i = 1; i < len; i++){
		if (s[sa[i]] != s[sa[i - 1]])m++;
		r[0][sa[i]] = m;
	}
	for (int i = 1; ++m < len; i *= 2){
		for (int j = 0; j < i; j++)
			tmp[j] = len - i + j;
		for (int j = i, k = 0; j < len; k++){
			if (sa[k] >= i)tmp[j++] = sa[k] - i;
		}
		memset(bucket, 0, sizeof(int)*m);
		for (int j = 0; j < len; j++)
			bucket[r1[j]]++;
		for (int j = 1; j < m; j++)
			bucket[j] += bucket[j - 1];
		for (int j = len - 1; j >= 0; j--)
			sa[--bucket[r1[tmp[j]]]] = tmp[j];
		r2[sa[0]] = m = 0;
		for (int j = 1; j < len; j++){
			if (r1[sa[j]] != r1[sa[j - 1]] || r1[sa[j] + i] != r1[sa[j - 1] + i])m++;
			r2[sa[j]] = m;
		}
		swap(r1, r2);
	}
	rk = r1;
	return len;
}
#else
int pool[3 * MAXN + 200], sa[3 * MAXN], rk[MAXN];
int tmp[MAXN];
void sort(int src[], int sa1[], int sa2[], int len, int m)
{
	memset(bucket, 0, sizeof(int)*m);
	for (int i = 0; i < len; i++)
		bucket[tmp[i] = src[sa1[i]]]++;
	for (int i = 1; i < m; i++)
		bucket[i] += bucket[i - 1];
	for (int i = len - 1; i >= 0; i--)
		sa2[--bucket[tmp[i]]] = sa1[i];
}
bool cmp(int r[], int r2[], int i, int j)
{
	if (r[i] != r[j])return r[i] < r[j];
	if (j % 3 == 1)return r2[i + 1] < r2[j + 1];
	if (r[i + 1] != r[j + 1])return r[i + 1] < r[j + 1];
	return r2[i + 2] < r2[j + 2];
}
#define F(pos) (pos % 3 == 1 ? pos / 3 : pos / 3 + t2)
#define G(pos) (pos < t2 ? pos * 3 + 1 : (pos - t2) * 3 + 2)
void merge(int r[], int sa[], int len, int m)
{
	int t1 = 0, t2 = (len + 1) / 3, t3 = len / 3, len2 = t2 + t3;
	int *r2 = r + len + 2, *sa2 = sa + len, *od = sa2 + len2;
	for (int i = 0, j = 0; i < len; i += 3, j += 2){
		sa[j] = i + 1;
		sa[j + 1] = i + 2;
	}
	r[len] = r[len + 1] = 0;
	sort(r + 2, sa, od, len2, m);
	sort(r + 1, od, sa, len2, m);
	sort(r, sa, od, len2, m);
	int k = 0;
	r2[F(od[0])] = 0;
	for (int i = 1; i < len2; i++){
		int pos1 = od[i], pos2 = od[i - 1];
		if (r[pos1] != r[pos2] || r[pos1 + 1] != r[pos2 + 1] || r[pos1 + 2] != r[pos2 + 2])k++;
		r2[F(od[i])] = k;
	}
	if (++k < len2)merge(r2, sa2, len2, k);
	else for (int i = 0; i < len2; i++)
		sa2[r2[i]] = i;
	if (len % 3 == 1)r2[t1++] = len - 1;
	for (int i = 0; i < len2; i++){
		if (sa2[i] < t2)r2[t1++] = sa2[i] * 3;
	}
	sort(r, r2, od, t1, m);
	for (int i = 0; i < len2; i++)
		r2[sa2[i] = G(sa2[i])] = i;
	int ii = 0, jj = 0, kk = 0;
	while (ii < t1 && jj < len2)
		sa[kk++] = cmp(r, r2, od[ii], sa2[jj]) ? od[ii++] : sa2[jj++];
	for (; ii < t1; ii++)sa[kk++] = od[ii];
	for (; jj < len2; jj++)sa[kk++] = sa2[jj];
}
int suffixArray()
{
	int len = 0;
	while (pool[len] = s[len], s[len++]);
	merge(pool, sa, len, CHAR);
	return len;
}
#endif
int h[MAXN];
int getH()
{
	int len = suffixArray() - 1;
#ifdef DC3
	for (int i = 0; i <= len; i++)
		rk[sa[i]] = i;
#endif
	for (int i = 0, k = 0; i < len; i++){
		int j = sa[rk[i] - 1];
		while (s[i + k] == s[j + k])k++;
		h[rk[i] - 1] = k;
		if (k)k--;
	}
	return len;
}
unsigned char lg2[MAXN * 2];
int st[20][MAXN];
void init(int len)
{
	for (int i = 0; i < len; i++)
		st[0][i] = h[i];
	for (unsigned char i = 0; (1 << i) <= len; i++)
		memset(&lg2[1 << i], i, 1 << i);
	for (int i = 1; i <= lg2[len]; i++){
		for (int k = len - (1 << i); k >= 0; k--)
			st[i][k] = min(st[i - 1][k], st[i - 1][k + (1 << (i - 1))]);
	}
}
inline int query(int i, int j){
	if (i > j)swap(i, j);
	int t = lg2[j - i];
	return min(st[t][i], st[t][j - (1 << t)]);
}
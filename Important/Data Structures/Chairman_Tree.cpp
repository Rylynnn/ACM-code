#include<cstdio>
#include<algorithm>
using namespace std;
#define MAXN 100001
struct Tree{
	int num, childL, childR;
}tree[MAXN * 20];
int root[MAXN];
int cnt, maxValue;
void init(int maxValue)
{
	root[0] = 0; cnt = 0;
	::maxValue = maxValue;
}
void insert(int& i, int l, int r, int value)
{
	tree[++cnt] = tree[i];
	tree[i = cnt].num++;
	int mid = (unsigned int)(l + r) >> 1;
	if (l == mid)return;
	if (value < mid)insert(tree[i].childL, l, mid, value);
	else insert(tree[i].childR, mid, r, value);
}
int query(int i, int j, int k, int l, int r)
{
	int mid = (unsigned int)(l + r) >> 1;
	if (l == mid)return l;
	int t = tree[tree[j].childL].num - tree[tree[i].childL].num;
	if (t >= k)return query(tree[i].childL, tree[j].childL, k, l, mid);
	else return query(tree[i].childR, tree[j].childR, k - t, mid, r);
}
int order[MAXN], a[MAXN];
int n, m;
int main()
{
	int test;
	scanf("%d", &test);
	for (int tt = 0; tt < test; tt++){
		scanf("%d%d", &n, &m);
		init(1 << 18);
		for (int i = 1; i <= n; i++)
			scanf("%d", &a[i]);
		memcpy(order, a + 1, sizeof(int)*n);
		sort(order, order + n);
		for (int i = 1; i <= n;i++){
			root[i] = root[i - 1];
			insert(root[i], 0, maxValue, lower_bound(order, order + n, a[i]) - order);
		}
		for (int t = 0; t < m; t++){
			int i, j, k;
			scanf("%d%d%d", &i, &j, &k);
			printf("%d\n", order[query(root[i - 1], root[j], k, 0, maxValue)]);
		}
	}
	return 0;
}
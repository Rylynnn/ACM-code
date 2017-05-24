#include<cstdio>
struct Int{
	int value;
	void add(int value){ this->value += value; }
	int sum(){ return value; }
};
template<typename T, int size> struct Tree{
	static int pos[3], cnt, bound;
	static Tree pool[size];
	T value;
	Tree *ch[2];
	void add(int value){
		Tree *rt = this;
		for (int i = bound; i;){
			i >>= 1; rt->value.add(value);
			int k = i & pos[2] ? 1 : 0;
			if (!rt->ch[k]){
				rt->ch[k] = ++cnt + pool;
				*rt->ch[k] = pool[0];
			}
			rt = rt->ch[k];
		}
	}
	int sum(){ return sum(0, bound); }
	int sum(int l, int r){
		if (l >= pos[0] && r <= pos[1])return value.sum();
		int mid = (l + r) >> 1, ret = 0;
		if (mid > pos[0] && ch[0])ret += ch[0]->sum(l, mid);
		if (mid < pos[1] && ch[1])ret += ch[1]->sum(mid, r);
		return ret;
	}
};
typedef Tree<Int, 100001 * 289> Tree1;
typedef Tree<Tree1, 100001 * 17> Tree2;
template<typename T, int s> int Tree<T, s>::pos[3];
template<typename T, int s> int Tree<T, s>::cnt;
template<typename T, int s> int Tree<T, s>::bound;
template<typename T, int s> Tree<T, s> Tree<T, s>::pool[s];
#include<algorithm>
#include<cstring>
#include<cmath>
using namespace std;
Tree2 rt;
int a[100001], b[100001], orderA[100001], orderB[100001];
struct Query{
	int pos, x, y, id;
	bool operator < (const Query& q)const{ return pos < q.pos; }
}q[200001];
int ans[200001];
int main()
{
	int n, m, test;
	scanf("%d", &test);
	while (test--){
		scanf("%d%d", &n, &m);
		for (int i = 0; i < n; i++)
			scanf("%d", &a[i]);
		for (int i = 0; i < n; i++)
			scanf("%d", &b[i]);
		memcpy(orderA, a, sizeof(int)*n);
		sort(orderA, orderA + n);
		memcpy(orderB, b, sizeof(int)*n);
		sort(orderB, orderB + n);
		for (int i = 0; i < n; i++){
			a[i] = lower_bound(orderA, orderA + n, a[i]) - orderA;
			b[i] = lower_bound(orderB, orderB + n, b[i]) - orderB;
		}
		for (int i = 0; i < m; i++){
			int l, r, x, y;
			scanf("%d%d%d%d", &l, &r, &x, &y);
			x = lower_bound(orderA, orderA + n, x) - orderA;
			y = lower_bound(orderB, orderB + n, y) - orderB;
			q[2 * i] = { l - 1, x, y, 2 * i };
			q[2 * i + 1] = { r, x, y, 2 * i + 1};
		}
		sort(q, q + 2 * m);
		rt = Tree2::pool[0];
		Tree1::cnt = Tree2::cnt = 0;
		Tree1::bound = Tree2::bound = 1 << (int)log2(2 * n - 1);
		Tree1::pos[1] = Tree2::pos[1] = n;
		for (int i = 0, j = 0; i < 2 * m; i++){
			for (; j < q[i].pos; j++){
				Tree1::pos[2] = a[j]; Tree2::pos[2] = b[j];
				rt.add(1);
			}
			Tree1::pos[0] = q[i].x; Tree2::pos[0] = q[i].y;
			ans[q[i].id] = rt.sum();
		}
		for (int i = 0; i < m; i++)
			printf("%d\n", ans[2 * i + 1] - ans[2 * i]);
	}
}
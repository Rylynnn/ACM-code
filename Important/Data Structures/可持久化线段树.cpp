struct Tree{
	long long value;
	int delta;
	int left, right;
}tree[325001];//空间大于2mlogn+n
int root[100001];
int cnt, op, n;
//初始化时cnt,op均为0;
//若初始化时a数组为0，则不用init，直接清空tree[0]。
void init(int a[], int n)
{
	tree[cnt].delta = 0;
	if (n == 1){ tree[cnt].value = *a; return; }
	int cur = cnt, mid = n >> 1;
	tree[cur].left = ++cnt;
	init(a, mid);
	tree[cur].right = ++cnt;
	init(a + mid, n - mid);
	tree[cur].value = tree[tree[cur].left].value + tree[tree[cur].right].value;
}
void pushDown(int i, int len)
{
	tree[++cnt] = tree[tree[i].left]; tree[i].left = cnt;
	tree[++cnt] = tree[tree[i].right]; tree[i].right = cnt;
	int t;
	if (t = tree[i].delta){
		tree[tree[i].left].value += (long long)t * (len / 2);
		tree[tree[i].left].delta += t;
		tree[tree[i].right].value += (long long)t * ((len + 1) / 2);
		tree[tree[i].right].delta += t;
		tree[i].delta = 0;
	}
}
int queryL, queryR, value;
void addInternal(int i, int l, int len)
{
	if (queryL <= l && queryR >= l + len){
		tree[i].value += (long long)value * len;
		tree[i].delta += value;
		return;
	}
	pushDown(i, len);
	int mid = l + len / 2;
	if (mid > queryL)addInternal(tree[i].left, l, len / 2);
	if (mid < queryR)addInternal(tree[i].right, mid, (len + 1) / 2);
	tree[i].value = tree[tree[i].left].value+ tree[tree[i].right].value;
}
//从0开始编号，不包括右端点，区间长度不能为0
void addValue(int l, int r, int num, int id){
	queryL = l; queryR = r; value = num;
	memcpy(&tree[++cnt], &tree[root[id]], sizeof(Tree));
	root[++op] = cnt;
	addInternal(cnt, 0, n);
}
long long queryInternal(int i, int l, int len)
{
	if (queryL <= l && queryR >= l + len)return tree[i].value;
	pushDown(i, len);
	int mid = l + len / 2;
	long long ret = 0;
	if (mid > queryL)ret += queryInternal(tree[i].left, l, len / 2);
	if (mid < queryR)ret += queryInternal(tree[i].right, mid, (len + 1) / 2);
	return ret;
}
//从0开始编号，不包括右端点，区间长度不能为0
inline long long query(int l, int r, int id){
	queryL = l; queryR = r;
	return queryInternal(root[id], 0, n);
}
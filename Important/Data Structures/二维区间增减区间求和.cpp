int n, m;
int tree[4][3002][3002];
int sum(int tree[][3002], int i, int j)
{
	int ret = 0;
	for (; i; i -= i&-i){
		for (int k = j; k; k -= k&-k)
			ret += tree[i][k];
	}
	return ret;
}
void add(int tree[][3002], int i, int j, int value)
{
	for (; i <= n; i += i&-i){
		for (int k = j; k <= m; k += k&-k)
			tree[i][k] += value;
	}
}
inline void F(int k, int i, int j, int value){
	int t = (k & 1 ? i : 1) * (k & 2 ? j : 1);
	add(tree[k], i, j, value*t);
}
void add(int i1, int j1, int i2, int j2, int value)
{
	for (int k = 0; k < 4; k++){
		F(k, i1, j1, value);
		F(k, i1, j2 + 1, -value);
		F(k, i2 + 1, j1, -value);
		F(k, i2 + 1, j2 + 1, value);
	}
}
inline int G(int k, int i, int j){
	int t = (k & 1 ? 1 : -(i + 1)) * (k & 2 ? 1 : -(j + 1));
	return sum(tree[k], i, j)*t;
}
int sum(int i1, int j1, int i2, int j2)
{
	int ret = 0;
	for (int k = 0; k < 4; k++){
		ret += G(k, i1 - 1, j1 - 1);
		ret -= G(k, i1 - 1, j2);
		ret -= G(k, i2, j1 - 1);
		ret += G(k, i2, j2);
	}
	return ret;
}

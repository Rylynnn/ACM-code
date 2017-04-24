#include<cstdio>
#include<algorithm>
using namespace std;
#define MAXN 300010
struct Tree{
	Tree *left, *right, *fa;
	int value, delta, max;
	bool reverse;
}pool[MAXN];
int num;//pool[0]代表空节点，判断指针是否为null要写rt!=root
inline void pushUp(Tree * rt){
	rt->max = max(rt->value, max(rt->left->max, rt->right->max));
}
inline void update(Tree *t, int delta){
	t->value += delta;
	t->delta += delta;
	t->max += delta;
}
inline void pushDown(Tree *rt){
	if (rt->delta){
		if (rt->left != pool)update(rt->left, rt->delta);
		if (rt->right != pool)update(rt->right, rt->delta);
		rt->delta = 0;
	}
	if (rt->reverse){
		rt->reverse = 0;
		rt->left->reverse ^= 1;
		rt->right->reverse ^= 1;
		swap(rt->left, rt->right);
	}
}
void rotate(Tree *t)
{
	Tree *p = t->fa, *q = p->fa;
	if (p->left == t){
		p->left = t->right;
		t->right->fa = p;
		t->right = p;
	}
	else{
		p->right = t->left;
		t->left->fa = p;
		t->left = p;
	}
	//此处和普通splay不同
	if (q->left == p)q->left = t;
	else if (q->right == p)q->right = t;
	p->fa = t; t->fa = q;
	pushUp(p);
}
inline bool check(Tree *t, Tree *f){ return f->left == t || f->right == t; }
void pushSign(Tree *t){
	Tree* f = t->fa;
	if (check(t, f))pushSign(f);
	pushDown(t);
}
void splay(Tree *t)
{
	pushSign(t);
	for (Tree* f; check(t, f = t->fa); rotate(t)){
		if (check(f, f->fa))
			rotate((f->left == t) == (f->fa->left == f) ? f : t);
	}
	pushUp(t);
}
Tree* access(Tree *rt)
{
	Tree *t;
	for (t = pool; rt != pool; rt = rt->fa){
		splay(rt);
		rt->right = t; t = rt;
		pushUp(t);
	}
	return t;
}
void setRoot(Tree *rt){ access(rt)->reverse ^= 1; }
Tree* findRoot(Tree *t)
{
	t = access(t); pushDown(t);
	for (; t->left != pool; pushDown(t))t = t->left;
	splay(t);
	return t;
}
bool connected(Tree *t1, Tree *t2)
{
	setRoot(t1);
	return findRoot(t2) == t1;
}
//有向树
void link(Tree *t1, Tree *t2){ splay(t1); t1->fa = t2; }
//无向树
void link(Tree *t1, Tree *t2)
{
	Tree *t = access(t1);
	t->reverse ^= 1; t->fa = t2;
}
//无向树，以t1为根拿出子树t2
void cut(Tree *t1, Tree *t2)
{
	setRoot(t1); splay(t2);
	if (t2->left){
		t2->left->fa = t2->fa;
		t2->left = t2->fa = pool;
		pushUp(t2);
	}
	else t2->fa = pool;
}
//有向树
void cut(Tree *t)
{
	splay(t);
	if (t->left){
		t->left->fa = t->fa;
		t->left = t->fa = pool;
	}
	else t->fa = pool;
}
void add(Tree *t1, Tree *t2, int delta)
{
	setRoot(t1);
	update(access(t2), delta);
}
int queryMax(Tree *t1, Tree *t2)
{
	setRoot(t1);
	return access(t2)->max;
}
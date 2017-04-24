#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
template<typename T>
struct Heap{
	T a[1000002];
	int n;
	Heap(){ clear(); }
	void down(int i){
		for (int j = i << 1; j <= n; i = j, j <<= 1){
			if (a[j + 1] < a[j])j++;
			if (a[j] < a[i])swap(a[i], a[j]);
			else break;
		}
	}
	void build(T src[], int n){
		memcpy(a + 1, src, n*sizeof(T));
		this->n = n;
		for (int i = n >> 1; i; i--){
			a[n + 1] = a[n];
			down(i);
		}
	}
	void push(T value){
		a[++n] = value;
		for (int j = n, i = j >> 1; i && a[j] < a[i]; j = i, i >>= 1)
			swap(a[i], a[j]);
	}
	void pop(){
		a[1] = a[n--]; 
		down(1);
	}
	T top()const{ return a[1]; }
	bool empty()const{ return n > 0; }
	void clear(){ n = 0; }
};
#include<cstdio>
#include<cstring>
using namespace std;
#define LETTER 26
struct Trie{
	int num, next;
}trie[1000000];
int cnt;
int pool[LETTER * 200000], poolEnd;
void init()
{
	cnt = 0;
	trie[0].num = 0;
	trie[0].next = -1;
	memset(pool, 0, sizeof(pool));
	poolEnd = 0;
}
inline int convert(char ch){ return ch - 'a'; }
inline char convert2(int value){ return value + 'a'; }
inline bool oneBranch(int value){ return value < LETTER; }
inline int child(int i, int ch){
	if (oneBranch(trie[i].next))return trie[i].next == ch ? i + 1 : 0;
	return pool[trie[i].next + ch];
}
void insert(char *s)
{
	int pos = 0, i;
	for (i = 0; s[i]; i++){
		int t = trie[pos].next;
		if (oneBranch(t)){
			if (t == convert(s[i]))pos++;
			else{
				trie[pos].next = (poolEnd += LETTER);
				if (t != -1)pool[trie[pos].next + t] = pos + 1;
				break;
			}
		}
		else if (pool[t + convert(s[i])])
			pos = pool[t + convert(s[i])];
		else break;
	}
	if (s[i]){
		pool[trie[pos].next + convert(s[i])] = ++cnt;
		for (i++; s[i]; i++, cnt++){
			trie[cnt].num = 0;
			trie[cnt].next = convert(s[i]);
		}
		trie[cnt].num = 1;
		trie[cnt].next = -1;
	}
	else trie[pos].num++;
}
int search(char* s)
{
	int pos = 0;
	for (int i = 0; s[i]; i++){
		pos = child(pos, convert(s[i]));
		if (!pos)return -1;
	}
	return trie[pos].num;
}
char temp[100];
void dfs(int i, int h)
{
	if (trie[i].num){
		temp[h] = 0;
		printf("%s %d\n", temp, trie[i].num);
	}
	int t = trie[i].next;
	if (oneBranch(t)){
		if (t == -1)return;
		temp[h] = convert2(t);
		dfs(i + 1, h + 1);
	}
	else for (int j = 0; j < LETTER; j++){
		if (pool[t + j]){
			temp[h] = convert2(j);
			dfs(pool[t + j], h + 1);
		}
	}
}
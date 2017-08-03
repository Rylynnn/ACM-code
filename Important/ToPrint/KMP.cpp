#include<cstdio>
#define MAXN 2000001
char s[MAXN];
int fail[MAXN];
void make_fail()
{
	for (int i = 1, j = 0; s[i]; i++){
		while (j && s[i] != s[j])j = fail[j - 1];
		if (s[i] == s[j])fail[i] = ++j;
		else fail[i] = 0;
	}
}
bool search(char *str)
{
	for (int i = 0, j = 0; str[i]; i++){
		while (j && str[i] != s[j])j = fail[j - 1];
		if (str[i] == s[j] && !s[++j])return true;
	}
	return false;
}
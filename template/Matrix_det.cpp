#include <iostream>
#include <cstring>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <queue>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <deque>
#include <bitset>
#include <algorithm>
using namespace std;
  
#define PI acos(-1.0)
#define LL long long
#define PII pair<int, int>
#define PLL pair<LL, LL>
#define mp make_pair
#define IN freopen("in.txt", "r", stdin)
#define OUT freopen("out.txt", "wb", stdout)
#define scan(x) scanf("%d", &x)
#define scan2(x, y) scanf("%d%d", &x, &y)
#define scan3(x, y, z) scanf("%d%d%d", &x, &y, &z)
#define sqr(x) (x) * (x)
#define pr(x) cout << #x << " = " << x << endl
#define lc o << 1
#define rc o << 1 | 1
#define pl() cout << endl
#define CLR(a, x) memset(a, x, sizeof(a))
#define FILL(a, n, x) for (int i = 0; i < n; i++) a[i] = x

LL det(LL a[][maxn], int n) {
    int i, j, k, r;
    LL res = 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            while (a[j][i]) {
                LL f = a[i][i] / a[j][i];
                for (int k = i; k < n; k++) a[i][k] -= f * a[j][k];
                for (int k = i; k < n; k++) swap(a[i][k], a[j][k]);
                res = -res;
            }
        }
        if (a[i][i] == 0) return 0;
        res *= a[i][i];
    }
    return res < 0 ? -res : res;
}


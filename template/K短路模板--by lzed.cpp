#include<iostream>
#include<cstring>
#include<stack>
#include<vector>
#include<set>
#include<map>
#include<algorithm>
#include<cmath>
#include<queue>
#include<sstream>
#include<iomanip>
#include<fstream>
#include<cstdio>
#include<cstdlib>
#include<climits>
#include<deque>
using namespace std;
  
#define INF 0x7fffffff
#define PI acos(-1.0)
#define LL long long
#define PII pair<int, int>
#define PLL pair<LL, LL>
#define mp make_pair
#define IN freopen("in.txt", "r", stdin)
#define OUT freopen("out.txt", "wb", stdout)
#define scan(x) scanf("%d", &x)
#define scan2(x, y) scanf("%d%d", &x, &y)
#define sqr(x) (x) * (x)

/*
 * 此做法只适用于无向图，若为有向图，考虑到h[]数组为某点到终点的最短距离，须建新图将所有有向边反向并进行SPFA操作
 */

const int maxn = 5001;
int N, R;
struct Edge {
    int from, to, dist;
    Edge(int u, int v, int w) : from(u), to(v), dist(w) {

    }
};

vector<Edge> edges;
vector<int> G[maxn];

void clr() {
    edges.clear();
    for (int i = 0; i < maxn; i++) G[i].clear();
}

void add(int u, int v, int w) {
    edges.push_back(Edge(u, v, w));
    int s = edges.size() - 1;
    G[u].push_back(s);
}


int h[maxn];
void spfa(int s) {
    int inq[maxn];
    memset(inq, 0, sizeof(inq));
    queue<int> q;
    for (int i = 1; i <= N; i++) h[i] = INF;
    h[s] = 0;
    inq[s] = 1;
    q.push(s);
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        inq[u] = 0;
        int s = G[u].size();
        for (int i = 0; i < s; i++) {
            Edge &e = edges[G[u][i]];
            if (h[u] < INF && h[e.to] > h[u] + e.dist) {
                h[e.to] = h[u] + e.dist;
                if (!inq[e.to]) {
                    inq[e.to] = 1;
                    q.push(e.to);
                }
            }
        }
    }
}

struct node {
    int _f, _g, _v;
    node (int x, int y, int z) : _f(x), _g(y), _v(z) {

    }
};

struct cmp {
    bool operator() (node x, node y) {
        return x._f > y._f;
    }
};

int K = 2;
int Astar(int s) {
    int cnt[maxn], shortest;
    memset(cnt, 0, sizeof(cnt));
    priority_queue<node, vector<node>, cmp> q;
    q.push(node(h[s], 0, s));
    while (!q.empty()) {
        node t = q.top();
        q.pop();
        int u = t._v;
        cnt[u]++;
        if (u == N && cnt[N] == K) return t._f;
        int s = G[u].size();
        for (int i = 0; i < s; i++) {
            Edge &e = edges[G[u][i]];
            int g = t._g + e.dist;
            int f = h[e.to] + g;
            q.push(node(f, g, e.to));
        }
    }
    return -1;
}

int main() {
    //IN;
    scan2(N, R);
    for (int i = 0; i < R; i++) {
        int x, y, z;
        scanf("%d%d%d", &x, &y, &z);
        add(x, y, z);
        add(y, x, z);
    }
    spfa(N);
    int res = Astar(1);
    cout << res << endl;
    return 0;
}



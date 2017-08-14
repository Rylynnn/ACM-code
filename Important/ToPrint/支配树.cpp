//Codechef May 15 GRAPHCNT
//一张n个点m条边的有向图，问有多少对（x,y）存在有一条从1到x的路径，有一条1到y的路径，他们俩不相交。
//求出支配树之后就是一个问树上多少点对lca是1
//如果对于r->w的任意一条路径中都存在一个点p，那么我们称点p为w的支配点（当然这也是r->w的必经点），注意r点不讨论支配点。下面用idom[u]表示离点u最近的支配点。 
//对于原图上除r外每一个点u，从idom[u]向u建一条边，最后我们可以得到一个以r为根的树。这个树我们就叫它“支配树”。
//
#include <iostream>
#include <cstdio>
using namespace std;
typedef long long lld;
const int MaxN = 100000 + 10, MaxE = (5 * 100000) * 2 + MaxN;
int head[MaxN], pre[MaxN], dom[MaxN], to[MaxE], nxt[MaxE], top;
void addedge(int *h,int fr,int tt)
{
    top ++;
    nxt[top] = h[fr];
    to[top] = tt;
    h[fr] = top;
}
int n, m;
void init()
{
    scanf("%d%d", &n, &m);
    int a, b;
    for(int i = 1; i <= m; ++i)
    {
        scanf("%d%d", &a, &b);
        addedge(head, a, b);
        addedge(pre, b, a);
    }
}
int bcj[MaxN], semi[MaxN], idom[MaxN], best[MaxN], dfn[MaxN], id[MaxN], fa[MaxN], dfs_clock;
int push(int v)
{
    if(v == bcj[v]) return v;
    int y = push(bcj[v]);
    if(dfn[semi[best[bcj[v]]]] < dfn[semi[best[v]]]) best[v] = best[bcj[v]];
    return bcj[v] = y;
}//带权并查集路径压缩
void dfs(int rt)
{
    dfn[rt] = ++dfs_clock;
    id[dfs_clock] = rt;
    for(int i = head[rt]; i; i = nxt[i])
        if(!dfn[to[i]])
        {
            dfs(to[i]);
            fa[to[i]] = rt;
        }

}//求出dfs序
void tarjan()
{
    for(int i = dfs_clock, u; i >= 2; --i)
    {//按dfs序从大到小计算
        u = id[i];
        for(int j = pre[u]; j; j = nxt[j])//semi
        {
            if(!dfn[to[j]]) continue;
            push(to[j]);
            if(dfn[semi[best[to[j]]]] < dfn[semi[u]]) semi[u] = semi[best[to[j]]];
        }
        addedge(dom, semi[u], u);
        bcj[u] = fa[u];u = id[i - 1];
        for(int j = dom[u]; j; j = nxt[j])//idom
        {
            push(to[j]);
            if(semi[best[to[j]]] == u) idom[to[j]] = u;
            else idom[to[j]] = best[to[j]];
        }
    }
    for(int i = 2, u; i <= dfs_clock; ++i)
    {
        u = id[i];
        if(idom[u] != semi[u]) idom[u] = idom[idom[u]];
    }
}
int sons[MaxN];
lld ans;
void calc_son()
{
    for(int i = dfs_clock, u; i >= 2; --i)
    {
        u = id[i];
        ++ sons[u];
        if(idom[u] != 1) sons[idom[u]] += sons[u];
        else ans -= ((lld)sons[u] * (lld)(sons[u] - 1)) / 2ll;
    }
}
void solve()
{
    for(int i = 1; i <= n; ++i) bcj[i] = semi[i] = best[i] = i;
    dfs_clock = 0;
    dfs(1);
    tarjan();
    ans = ((lld)dfs_clock * (lld)(dfs_clock - 1)) / 2ll;
    calc_son();
    cout << ans << endl;
}
int main()
{
    init();
    solve();
    return 0;
}

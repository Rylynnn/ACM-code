struct MDST
{
    int n;
    int w[maxn][maxn];
    int vis[maxn];
    int res;
    int removed[maxn];
    int cid[maxn];
    int pre[maxn];
    int in[maxn];
    int max_cid;

    void init(int n)
    {
        this->n = n;
        memset(w,INF,sizeof(w));
    }

    void addEdge(int u,int v,int cost)
    {
        w[u][v] = min(w[u][v],cost);
    }

    int dfs(int u)
    {
        vis[u] = 1;
        int cnt = 1;
        for(int i = 1;i <= n;++i)
            if(!vis[i] && w[u][i] < INF) cnt += dfs(i);
        return cnt;
    }

    bool cycle(int u)
    {
        max_cid++;
        int v = u;
        while(cid[v] != max_cid) { cid[v] = max_cid; v = pre[v]; }
        return v == u;
    }

    void update(int u)
    {
        in[u] = INF;
        for(int i = 1;i <= n;++i)
            if(!removed[i] && w[i][u] < in[u])
            {
                in[u]= w[i][u];
                pre[u] = i;
            }
    }


    bool getRes(int s)
    {
        memset(vis,0,sizeof(vis));
        if(dfs(s) != n) return false;

        memset(removed,0,sizeof(removed));
        memset(cid,0,sizeof(cid));
        for(int i = 1;i <= n;++i) update(i);
        pre[s] = s,in[s] = 0;
        res = max_cid = 0;
        while(1)
        {
            bool have_cycle = false;
            for(int u = 1;u <= n;++u)
                if(u != s && !removed[u] && cycle(u))
                {
                    have_cycle = true;
                    int v = u;
                    do
                    {
                        if(v != u) removed[v] = 1;
                        res += in[v];

                        for(int i = 1;i <= n;++i)
                            if(cid[i] != cid[u] && !removed[i])
                            {
                                if(w[i][v] < INF) w[i][u] =
                                    min(w[i][u],w[i][v]-in[v]);
                                w[u][i] = min(w[u][i],w[v][i]);
                                if(pre[i] == v) pre[i] = u;
                            }
                        v = pre[v];
                    }while(v != u);
                    update(u);
                    break;
                }
                if(!have_cycle) break;
        }
        for(int i = 1;i <= n;++i)
            if(!removed[i]) res += in[i];
        return true;
    }
};


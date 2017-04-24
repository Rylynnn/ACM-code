void cal_sz(int v,int p)
{
    sz[v]=1;
    for(auto x : g[v])
    {
        if(x!=p)
        {
            cal_sz(x,v);
            sz[v]+=sz[x];
        }
    }
}

void add(int v,int p,int x)
{
    add_szsz(cnt[tr[v].col],-1);
    cnt[tr[v].col]+=x;
    add_szsz(cnt[tr[v].col],1);
    for(auto u : g[v])
    {
        if(u!=p && !big[u])
            add(u,v,x);
    }
}

void DFS(int v,int p,bool k)
{
    int bn=-1;
    for(auto u : g[v])
    {
        if(u!=p)
        {
            if(bn==-1 || sz[u]>sz[bn])bn=u;
        }
    }
    for(auto u : g[v])
    {
        if(u!=p && u!=bn)
            DFS(u,v,false);
    }
    if(bn!=-1)
    {
        DFS(bn,v,true);
        big[bn]=true;
    }
    add(v,p,1);
    //now cnt is the total of subtree of node v(include v),you can query now
    if(bn!=-1)big[bn]=false;
    if(!k)add(v,p,-1);
}

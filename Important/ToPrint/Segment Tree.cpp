struct segtree
{
    struct node
    {
        int ls,rs,lg,rg,mid;
        LL mi,ma,tag;
        void cle(){ls=rs=lg=rg=mid=mi=ma=tag=0;}
        bool one(){return lg==rg;}
        bool judge(int x,int y){return lg==x && rg==y;}
    }tr[200005];
    int k;

    void pushup(int x)
    {
        int ls=tr[x].ls,rs=tr[x].rs;
        if(ls+rs<=0)return;
        tr[x].mi=min(tr[ls].mi,tr[rs].mi);
        tr[x].ma=max(tr[ls].ma,tr[rs].ma);
    }
    
    void pushdown(int x)
    {
        int ls=tr[x].ls,rs=tr[x].rs;
        if(ls+rs<=0 || tr[x].tag==0)return;
        LL tag=tr[x].tag;
        tr[x].tag=0;
        tr[ls].tag+=tag;tr[rs].tag+=tag;
        tr[ls].mi+=tag;tr[ls].ma+=tag;
        tr[rs].mi+=tag;tr[rs].ma+=tag;
        pushup(x);
    }

    void mt(int x,int y)
    {
        tr[k].cle();
        tr[k].lg=x;tr[k].rg=y;
        if(x==y)return;
        int mid=(x+y)>>1,t=k;
        tr[k].mid=mid;
        k++;tr[t].ls=k;mt(x,mid);
        k++;tr[t].rs=k;mt(mid+1,y);
        pushup(t);
    }
    
    void init(int n){k=1;mt(0,n);}

    void add(int now,int x,int y,LL nu)
    {
        if(tr[now].judge(x,y))
        {
            tr[now].tag+=nu;
            tr[now].ma+=nu;
            tr[now].mi+=nu;
            return;
        }
        pushdown(now);
        int mid=tr[now].mid;
        if(x<=mid)add(tr[now].ls,x,min(mid,y),nu);
        if(y>mid)add(tr[now].rs,max(x,mid+1),y,nu);
        pushup(now);
    }

    pll query(int now,int x,int y)
    {
        if(tr[now].judge(x,y))return mp(tr[now].mi,tr[now].ma);
        int mid=tr[now].mid;
        pushdown(now);
        if(y<=mid)return query(tr[now].ls,x,y);
        if(x>mid)return query(tr[now].rs,x,y);
        pll u=query(tr[now].ls,x,mid);
        pll v=query(tr[now].rs,mid+1,y);
        return mp(min(u.fir,v.fir),max(u.sec,v.sec));
    }
}segt;

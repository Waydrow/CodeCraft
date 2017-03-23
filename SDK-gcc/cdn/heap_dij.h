//..xt:bfs���зֲ�
bool bfs()
{
    static int Q[MAXN];
    memset(sp,-1,sizeof(sp));
    sp[S]=0; Q[0]=S;
    edge *it;
    for (int h=0,t=1,u,v;h<t;++h)
    {
        for (u=Q[h],it=V[u];it;it=it->next)
        {
            v=it->t;
            if (sp[v]==-1&&it->c>0)
            {
                sp[v]=sp[u]+1; Q[t++]=v;
            }
        }
    }
    return sp[T]!=-1;
}
//..xt:dfs������·
int dfs(int u,int low)
{
    if (u==T) return low;
    int ret=0,tmp,v;
    for (edge *&it=cur[u];(it)&&ret<low;it=it->next)
    {
        v=it->t;
        if (sp[v]==sp[u]+1&&it->c>0)
        {
            if ((tmp=dfs(v,min(low-ret,it->c)))>0)
            {
                ret+=tmp; it->c-=tmp; it->op->c+=tmp;
            }
        }
    }
    if (!ret) sp[u]=-1; return ret;
}

//..xt:ÿ���ҵ�
//..xt������ͼ�Ժ󷵻���С���������
int maxcostflow() {
    int tmp;
    while(bfs()) {
        memcpy(cur,V,sizeof(V));
        while((tmp=dfs(S,INF))>0);
    }
    return 0;
}

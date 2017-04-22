//..xt:如果存在从源点到汇点大路径则返回true，否则返回false
bool SPFA() {
    int u,v;
    for(u=S; u<=T; u++) {
        sp[u]=INF;
    }
    int q[MAXN];
    int st=0,ed=0;
    Prev[S]=-1;
    q[ed++]=S;
    sp[S]=0;
    vis[S]=1;
    int sum=0,num=1;
    double ave=0;
    while(st!=ed) {
        u=q[st++];
        if(st>=MAXN)st=0;
        while(sp[u]>ave){
            q[ed++]=u;
            if(ed>=MAXN)ed=0;
            u=q[st++];
            if(st>=MAXN)st=0;
        }
        num--;
        sum-=sp[u];
        if(num)ave=(double)sum/num;
        else ave=0;
        vis[u]=0;
        for(edge *k=V[u]; k; k=k->next) {
            v=k->t;
            if(k->c>0&&sp[u]+k->v<sp[v]) {
                sp[v]=sp[u]+k->v;
                Prev[v]=u;
                path[v]=k;
                if(vis[v]==0) {
                    if(num==0||sp[v]>sp[q[st]]){
                        q[ed++]=v;
                        if(ed>=MAXN)ed=0;
                    }
                    else{
                        st--;
                        if(st<0)st=MAXN-1;
                        q[st]=v;
                    }
                    num++;
                    sum+=sp[v];
                    ave=(double)sum/num;
                    vis[v]=1;
                }
            }
        }
    }
    return sp[T]!=INF;
}
//..xt:每次找到最短路进行增广，返回该增光路cost
int argument() {
    int i,cost=INF,flow=0;
    edge *e;
    for(i=T; Prev[i]!=-1; i=Prev[i]) {
        e=path[i];
        if(e->c<cost)cost=e->c;
    }
    for(int i=T; Prev[i]!=-1; i=Prev[i]) {
        e=path[i];
        e->c-=cost;
        e->op->c+=cost;
        flow+=e->v*cost;
    }
    return flow;
}
//..xt：建好图以后返回最小费用最大流
int maxcostflow() {
    int Flow=0;
    while(SPFA()) {
        Flow+=argument();
    }
    return Flow;
}

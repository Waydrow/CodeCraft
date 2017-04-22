vector<int> curRoad;//寻找路径dfs时记录
vector<vector<int> >relRoad;//输出答案，记录所有路径
vector<int>mustChoose;//记录一定需要选择的点的编号
map<unsigned long long,int>mp;//记录已经计算过的答案
// 最终需要输出到文件的内容
char topo_file[20000];
//..xt:费用流用到全局变量以及结构体
bool vis[MAXN];
struct edge {
    edge *next,*op;
    int t,c,v,bf;
} ES[MAXM],*V[MAXN],*cur[MAXN];
int N,M,S,T,EC=-1;
int demond[MAXN],sp[MAXN],Prev[MAXN];
edge *path[MAXN];
//记录链路信息
int links[50500][4];
int consumptionNodes[505][3];
int myNodes[1000];
int maxcostflow();

//堆优化dij
struct st{
    int id,dis;
};
bool operator <(const st &a,const st &b){
        return a.dis>b.dis;
    }
bool h_dij(){
    priority_queue<st>q;
    int u,v;
    for(u=S; u<=T; u++) {
        sp[u]=INF;
    }
    Prev[S]=-1;
    q.push({S,0});
    while(!q.empty()){
        st tmp=q.top();
        q.pop();
        u=tmp.id;
        int dis=tmp.dis;
        //printf("out id=%d,dis=%d\n",u,dis);
        if(dis>sp[u])continue;
        sp[u]=dis;
        for(edge *k=V[u]; k; k=k->next){
            v=k->t;
            if(k->c>0&&sp[u]+k->v<sp[v]) {
                sp[v]=sp[u]+k->v;
                Prev[v]=u;
                path[v]=k;
                //printf("in id=%d,dis=%d\n",v,sp[v]);
                q.push({v,sp[v]});
            }
        }
    }
    return sp[T]!=INF;
}
//..xt:将bitset进行hash存取
unsigned long long getHash(bitset<BITSIZE>&gene) {
    char hashStr[2000];
    int tmpChar=0;
    int pos=0;
    for(int i=0; i<nodesNum; i++) {
        tmpChar<<=1;
        tmpChar|=gene[i];
        if(i%4==0) {
            hashStr[pos++]=(char)(tmpChar+'a'-1);
            tmpChar=0;
        }
    }
    if(tmpChar)hashStr[pos++]=(char)(tmpChar+'a'-1);
    unsigned long long seed = 131; // 31 131 1313 13131 131313 etc..
    unsigned long long hash = 0;
    char *str=hashStr;
    while (*str) {
        hash = hash * seed + (*str++);
    }
    return hash;
}
//..xt:增加边的函数
void addedge(int a,int b,int v,int c=INF) {
    edge e1= {V[a],0,b,c,v,c},e2= {V[b],0,a,0,-v,0};
    ES[++EC]=e1;
    V[a]=&ES[EC];
    ES[++EC]=e2;
    V[b]=&ES[EC];
    V[a]->op=V[b];
    V[b]->op=V[a];
}
//..xt:读取数据
void readData(char *topo[],int nodesNum,int linkNum,int clientNum) {
    int lineAdd=4,numAdd=1;
    int st,ed,bandW,cost;
    for(int i=0; i<linkNum; i++) {
        sscanf(topo[i+lineAdd],"%d%d%d%d",&st,&ed,&bandW,&cost);
        links[i][0]=st;
        links[i][1]=ed;
        links[i][2]=bandW;
        links[i][3]=cost;
    }
    lineAdd=linkNum+5;
    numAdd+=nodesNum;
    int n_client,n_node,need;
    for(int i=0; i<clientNum; i++) {
        sscanf(topo[i+lineAdd],"%d%d%d",&n_client,&n_node,&need);
        myNodes[n_client]=n_node;
        consumptionNodes[i][0]=n_client;
        consumptionNodes[i][1]=n_node;
        consumptionNodes[i][2]=need;
    }
}
//..xt:建图函数
//..xt:网络节点编号：1...nodesNum
//..xt:消费节点编号：nodesNum+1...nodesNum+clientNum
//..xt:超级源点编号为0,超级汇点编号为nodesNum+clientNum+1
void buildBasicGraph(int nodesNum,int linkNum,int clientNum) {
    if(EC>=0) {
        for(int i=0; i<=EC; i++) {
            ES[i].c=ES[i].bf;
        }
        return;
    }
    S=0;
    T=nodesNum+clientNum+1;
    memset(V,0,sizeof(V));
    EC=-1;
    int numAdd=1;
    int st,ed,bandW,cost;
    for(int i=0; i<nodesNum; i++) {
        addedge(S,i+1,0,0);
    }
    for(int i=0; i<linkNum; i++) {
        st=links[i][0];
        ed=links[i][1];
        bandW=links[i][2];
        cost=links[i][3];
        st+=numAdd;
        ed+=numAdd;
        addedge(st,ed,cost,bandW);
        addedge(ed,st,cost,bandW);

    }
    numAdd+=nodesNum;
    int n_client,n_node,need;
    for(int i=0; i<clientNum; i++) {
        n_client=consumptionNodes[i][0];
        n_node=consumptionNodes[i][1];
        need=consumptionNodes[i][2];
        n_node++;
        n_client+=numAdd;
        addedge(n_node,n_client,0);
        addedge(n_client,T,0,need);
    }
}
//..xt:改变某条链路的容量
void setCapacity(int st,int ed,int c) {
    for(edge *k=V[st]; k; k=k->next) {
        if(k->t==ed) {
            k->c=c;
            k->op->c=0;
            break;
        }
    }
}
//..xt:检测是不是所有消费节点都得到满足,返回未满足的个数,并在imNotOk中做标记
bool imNotOk[MAXN];
int checkSatisfy(int nodesNum,int clientNum) {
    memset(imNotOk,0,sizeof(imNotOk));
    int numAdd=1+nodesNum;
    int rel=0;
    for(int i=0; i<clientNum; i++) {
        int u=i+numAdd;
        for(edge *k=V[u]; k; k=k->next) {
            if((k->t==T)&&(k->c!=0)){
                rel++;
                imNotOk[u]=true;
            }
        }
    }
    return rel;
}
//..xt:将流需求得不到满足的节点的流操作的费用进行撤销，算法为随机撤销，不考虑费用
int extraCost;
int sumCost;
int getExtraCost(int pos,int flow){
    if(pos>nodesNum){
        if(imNotOk[pos]){
            extraCost+=flow*sumCost;
        }
        return flow;
    }
    int totalFlow=0;
    for(edge *i=V[pos];i;i=i->next){
        if(i->c < i->bf){
            sumCost+=i->v;
            int tmpFlow=getExtraCost(i->t,min(flow,i->bf-i->c));
            sumCost-=i->v;
            totalFlow+=tmpFlow;
            i->c+=tmpFlow;
        }
    }
    return totalFlow;
}
//..xt:输出结果之前得到真实的relDNA(而不是我们认为这里存在服务器部署)
bitset<BITSIZE> getBetter() {
    int numAdd=1+nodesNum;
    bitset<BITSIZE>rel;
    rel.reset();
    for(int i=0; i<clientNum; i++) {
        int u=i+numAdd;
        for(edge *k=V[u]; k; k=k->next) {
            if((k->t==T)&&(k->c!=0)) rel[myNodes[i]]=1;
        }
    }
    return rel;
}
//..xt:给出某个体的评估值(严格使用最小费用最大流评估)
long long mynum=0;
long long mytime=0;
int calCost(bitset<BITSIZE>gene,int mustCal,bool deleteExtraCost) {
    while(1)puts("ggg");
    unsigned long long relHash=getHash(gene);
    if((!mustCal)&&mp.find(relHash)!=mp.end())return mp[relHash];
    buildBasicGraph(nodesNum,linkNum,clientNum);
    int numAdd=1;

    for(int i=0; i<nodesNum; i++) {
        if(gene[i]) {
            setCapacity(S,i+numAdd,INF);
        } else {
            setCapacity(S,i+numAdd,0);
        }
    }
    int sst=clock();
    int costRel=maxcostflow();
    cout <<"Normal Cost: "<<costRel<<endl;
    int eed=clock();
    int addit=checkSatisfy(nodesNum,clientNum);

    mynum++;
    mytime+=(eed-sst);
    extraCost=0;
    sumCost=0;
    if(deleteExtraCost){
        for(int i=1;i<=nodesNum;i++){
            getExtraCost(i,INF);
        }
    }
    return mp[relHash]=costRel+(gene.count()+addit)*deployCost-extraCost;
}
//..xt:递归寻找路径（用于打印结果）
int findRoad(int pos,int flow) {
    if(flow==0)return 0;
    if(pos>nodesNum) {
        curRoad.push_back(pos-1-nodesNum);
        curRoad.push_back(flow);
        relRoad.push_back(curRoad);
        curRoad.pop_back();
        curRoad.pop_back();
        return flow;
    }
    curRoad.push_back(pos-1);
    int sumFlow=0;
    for(edge *k=V[pos]; k; k=k->next) {
        if(k->c < k->bf) {
            int tmpFlow=findRoad(k->t,min(flow,k->bf - k->c));
            k->c+=tmpFlow;
            flow-=tmpFlow;
            sumFlow+=tmpFlow;
        }
    }
    curRoad.pop_back();
    return sumFlow;
}
//..xt:打印结果函数
void printRel(bitset<BITSIZE>rel) {
    for(int i=0; i<nodesNum; i++) {
        if(rel[i]) {
            findRoad(i+1,INF);
        }
    }
    int siz=relRoad.size();
    sprintf(topo_file, "%d\n\n", siz);
    relStr+=string(topo_file);
    for(int i=0; i<siz; i++) {
        int mm=relRoad[i].size();
        for(int j=0; j<mm; j++) {
            sprintf(topo_file, "%d ", relRoad[i][j]);
            relStr+=string(topo_file);
        }
        if(i<siz-1) {
            sprintf(topo_file, "\n");
            relStr+=string(topo_file);
        }
    }
}
//..xt：打印基因用于调试
void printGene(bitset<BITSIZE>gene) {
    for(int i=0; i<nodesNum; i++) {
        cout << gene[i];
    }
    cout << endl;
}
//..xt:贪心找出一定需要选择的直接连接的网络节点
void getMustChoose(){
    vector<pair<int ,int> > mp[1050];
    int num[1050];
    memset(num,0,sizeof(num));
    for(int i=0;i<linkNum;i++){
        int a=links[i][0];
        int b=links[i][1];
        int c=links[i][2];
        int v=links[i][3];
        num[a]+=c;
        num[b]+=c;
        mp[a].push_back(make_pair(v,c));
        mp[b].push_back(make_pair(v,c));
    }
    for(int i=0;i<clientNum;i++){
        int a=consumptionNodes[i][1];
        int c=consumptionNodes[i][2];
        if(num[a]<c){
            mustChoose.push_back(a);
            continue;
        }
        if(mp[a].size()){
            sort(mp[a].begin(),mp[a].end());
            int sum=0,siz=mp[a].size(),fee=0;
            for(int j=0;j<siz;j++){
                fee+=(mp[a][j].first)*min(c-sum,mp[a][j].second);
                sum+=min(c-sum,mp[a][j].second);
                if(sum>=c)break;
            }
            if(fee>=deployCost)mustChoose.push_back(a);
        }
    }
}


vector<int> curRoad;//Ѱ��·��dfsʱ��¼
vector<vector<int> >relRoad;//�����𰸣���¼����·��
vector<int>mustChoose;//��¼һ����Ҫѡ���ĵ��ı���
map<unsigned long long,int>mp;//��¼�Ѿ��������Ĵ���
// ������Ҫ�������ļ�������
char topo_file[20000];
//..xt:�������õ�ȫ�ֱ����Լ��ṹ��
bool vis[MAXN];
struct edge {
    edge *next,*op;
    int t,c,v,bf;
} ES[MAXM],*V[MAXN],*cur[MAXN];
int N,M,S,T,EC=-1;
int demond[MAXN],sp[MAXN],Prev[MAXN];
edge *path[MAXN];
//��¼��·��Ϣ
int links[50500][4];
int consumptionNodes[505][3];
int myNodes[1000];
int maxcostflow();

//���Ż�dij
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
//..xt:��bitset����hash��ȡ
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
//..xt:���ӱߵĺ���
void addedge(int a,int b,int v,int c=INF) {
    if(a==0&&gene[b-1]){
        out<<"a "<<a+1<<" "<<b+1<<" 0 "<<INF<<" "<<v<<endl;
    }
    else if(a!=0&&b!=nodesNum+clientNum+1){
        out<<"a "<<a+1<<" "<<b+1<<" 0 "<<c<<" "<<v<<endl;
    }
    edge e1= {V[a],0,b,c,v,c},e2= {V[b],0,a,0,-v,0};
    ES[++EC]=e1;
    V[a]=&ES[EC];
    ES[++EC]=e2;
    V[b]=&ES[EC];
    V[a]->op=V[b];
    V[b]->op=V[a];
}
//..xt:��ȡ����
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
//..xt:��ͼ����
//..xt:�����ڵ����ţ�1...nodesNum
//..xt:���ѽڵ����ţ�nodesNum+1...nodesNum+clientNum
//..xt:����Դ������Ϊ0,������������ΪnodesNum+clientNum+1
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
//..xt:�ı�ĳ����·������
void setCapacity(int st,int ed,int c) {
    for(edge *k=V[st]; k; k=k->next) {
        if(k->t==ed) {
            k->c=c;
            k->op->c=0;
            break;
        }
    }
}
//..xt:�����ǲ����������ѽڵ㶼�õ�����,����δ�����ĸ���,����imNotOk��������
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
//..xt:���������ò��������Ľڵ����������ķ��ý��г������㷨Ϊ���������������Ƿ���
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
//..xt:��������֮ǰ�õ���ʵ��relDNA(������������Ϊ�������ڷ���������)
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
//..xt:����ĳ����������ֵ(�ϸ�ʹ����С��������������)
long long mynum=0;
long long mytime=0;
int calCost(bitset<BITSIZE>gene,int mustCal,bool deleteExtraCost) {
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
    return costRel;
    /*
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
    int ddd=costRel+(gene.count()+addit)*deployCost-extraCost;
    printf("%d\n",ddd-(gene.count()+addit)*deployCost);
    return mp[relHash]=ddd;
    */
}
//..xt:�ݹ�Ѱ��·�������ڴ�ӡ������
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

//..xt����ӡ�������ڵ���
void printGene(bitset<BITSIZE>gene) {
    for(int i=0; i<nodesNum; i++) {
        cout << gene[i];
    }
    cout << endl;
}
//..xt:̰���ҳ�һ����Ҫѡ����ֱ�����ӵ������ڵ�
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

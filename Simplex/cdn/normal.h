vector<int> curRoad;//Ѱ��·��dfsʱ��¼
vector<vector<int> >relRoad;//�����𰸣���¼����·��
vector<int>mustChoose;//��¼һ����Ҫѡ���ĵ��ı���
set<int>notChoose;
map<unsigned long long,int>mp;//��¼�Ѿ��������Ĵ���


bool vis[MAXN];

int N,M,S,T,EC=-1;
int demond[MAXN],sp[MAXN],Prev[MAXN];
edge *path[MAXN];

int links[50500][4];
int consumptionNodes[505][3];
int myNodes[BITSIZE];
int serverInfo[15][2];
int nodesMoney[BITSIZE];
vector<int> NodesNeighbor[BITSIZE];
vector<int> roadCost[BITSIZE];
vector<int> roadBanw[BITSIZE];

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
    char hashStr[BITSIZE];
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
void addedge(int a,int b,int v,int c) {
    myEdgeNum[a]++;
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
    serverLevel=0;
    dataProperty datap=num_line;
    int tmpa,tmpb,tmpc,tmpd;
    int linksI=0,consumptionNodesI=0;
    for(int i=0;;i++){
        if(consumptionNodesI==clientNum){
            break;
        }
        if(topo[i][0]=='\r'||topo[i][0]=='\n'){
            datap = dataProperty(datap + 1);
            continue;
        }
        switch(datap){
        case num_line:
            continue;
        case servers_line:
            serverLevel++;
            sscanf(topo[i],"%d%d%d",&tmpa,&tmpb,&tmpc);
            serverInfo[tmpa][0]=tmpb;
            serverInfo[tmpa][1]=tmpc;
            break;
        case nodes_line:
            sscanf(topo[i],"%d%d",&tmpa,&tmpb);
            nodesMoney[tmpa]=tmpb;
            break;
        case links_line:
            sscanf(topo[i],"%d%d%d%d",&tmpa,&tmpb,&tmpc,&tmpd);
            links[linksI][0]=tmpa;
            links[linksI][1]=tmpb;
            links[linksI][2]=tmpc;
            links[linksI][3]=tmpd;
            NodesNeighbor[tmpa].push_back(tmpb);
            roadBanw[tmpa].push_back(tmpc);
            roadCost[tmpa].push_back(tmpd);
            NodesNeighbor[tmpb].push_back(tmpa);
            roadBanw[tmpb].push_back(tmpc);
            roadCost[tmpb].push_back(tmpd);
            linksI++;
            break;
        case consumtionNodes_line:
            sscanf(topo[i],"%d%d%d",&tmpa,&tmpb,&tmpc);
            consumptionNodes[consumptionNodesI][0]=tmpa;
            consumptionNodes[consumptionNodesI][1]=tmpb;
            consumptionNodes[consumptionNodesI][2]=tmpc;
            myNodes[tmpa]=tmpb;
            consumptionNodesI++;
        }
    }
}

void buildBasicGraph() {
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
        addedge(n_node,n_client,0,need);
        addedge(n_client,T,0,need);
    }
}

void setCapacity(int st,int ed,int c) {
    for(edge *k=V[st]; k; k=k->next) {
        if(k->t==ed) {
            k->c=c;
            k->op->c=0;
            break;
        }
    }
}

bool imNotOk[MAXN];
int checkSatisfy(bitset<BITSIZE>gene,int nodesNum,int clientNum) {
    vector<pair<int,int> >rrr;
    memset(imNotOk,0,sizeof(imNotOk));
    int numAdd=1+nodesNum;
    int rel=0;
    for(int i=0; i<clientNum; i++) {
        int u=i+numAdd;
        for(edge *k=V[u]; k; k=k->next) {
            if((k->t==T)&&(k->c!=0)){
                if(gene[ myNodes[i] ]){
                    while(1){
                        puts("adderror1");
                    }
                }
                for(int w=0;w<serverLevel;w++){
                    if(serverInfo[w][0]>=(k->c)){
                        rel+=serverInfo[w][1];
                        rel+=nodesMoney[ myNodes[i] ];
                        break;
                    }
                    if(w==serverLevel-1){
                        while(1){
                            puts("adderror2");
                        }
                    }
                }
                imNotOk[u]=true;
            }
        }
    }
    return rel;
}
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
long long mynum=0;
long long mytime=0;
int calCost(bitset<BITSIZE>gene,int mustCal,bool deleteExtraCost) {
    unsigned long long relHash=getHash(gene);
    if((!mustCal)&&mp.find(relHash)!=mp.end())return mp[relHash];
    buildBasicGraph();
    int numAdd=1;
    int rel=0;
    for(int i=0; i<nodesNum; i++) {
        if(gene[i]) {
            setCapacity(S,i+numAdd,serverInfo[serverLevel-1][0]);
        } else {
            setCapacity(S,i+numAdd,0);
        }
    }
    int sst=clock();
    int costRel=maxcostflow();
    rel+=costRel;
    int myflow[BITSIZE];
    memset(myflow,0,sizeof(myflow));
    for(int i=1;i<=nodesNum;i++){
        for(edge *k=V[i];k;k=k->next){
            if(!(k->bf))continue;
            int tt=k->t;
            myflow[i-1]+=(k->bf-k->c);
            if(tt-1<nodesNum){
                myflow[tt-1]-=(k->bf-k->c);
            }
        }
    }
    for(int i=0;i<nodesNum;i++){
        while(myflow[i]<0||myflow[i]>serverInfo[serverLevel-1][0]){
            puts("flowadderror");
        }
        if(gene[i]){
            for(int w=0;w<serverLevel;w++){
                if(serverInfo[w][0]>=myflow[i]){
                    rel+=serverInfo[w][1];
                    break;
                }
                while(w==serverLevel-1){
                    puts("adderror3");
                }
            }
        }

    }
    for(int i=0;i<nodesNum;i++){
        if(gene[i]){
            rel+=nodesMoney[i];
        }
    }
    int eed=clock();
    int addit=checkSatisfy(gene,nodesNum,clientNum);
    rel+=addit;
    mynum++;
    mytime+=(eed-sst);
    extraCost=0;
    sumCost=0;
    /*
    if(deleteExtraCost){
        for(int i=1;i<=nodesNum;i++){
            getExtraCost(i,INF);
        }
    }
    */
    //printf("rrr: %d\n",rrr);
    return mp[relHash]=rel-extraCost;
}
//..xt:�ݹ�Ѱ��·�������ڴ�ӡ������
int findRoad(int pos,int flow,int id) {
    if(flow==0)return 0;
    if(pos>nodesNum) {
        curRoad.push_back(pos-1-nodesNum);
        curRoad.push_back(flow);
        curRoad.push_back(id);
        relRoad.push_back(curRoad);
        curRoad.pop_back();
        curRoad.pop_back();
        curRoad.pop_back();
        return flow;
    }
    curRoad.push_back(pos-1);
    int sumFlow=0;
    for(edge *k=V[pos]; k; k=k->next) {
        if(k->c < k->bf) {
            int tmpFlow=findRoad(k->t,min(flow,k->bf - k->c),id);
            k->c+=tmpFlow;
            flow-=tmpFlow;
            sumFlow+=tmpFlow;
        }
    }
    curRoad.pop_back();
    return sumFlow;
}
//..xt:��ӡ��������
void printRel(bitset<BITSIZE>rel) {
    int myflow[BITSIZE];
    memset(myflow,0,sizeof(myflow));
    for(int i=1;i<=nodesNum;i++){
        for(edge *k=V[i];k;k=k->next){
            if(!(k->bf))continue;
            int tt=k->t;
            myflow[i-1]+=(k->bf-k->c);
            if(tt-1<nodesNum){
                myflow[tt-1]-=(k->bf-k->c);
            }
        }
    }
    for(int i=0; i<nodesNum; i++) {
        if(rel[i]) {
            int flowz=0,id=0;
            for(int w=0;w<serverLevel;w++){
                if(serverInfo[w][0]>=myflow[i]){
                    flowz=serverInfo[w][0];
                    id=w;
                    break;
                }
            }
            findRoad(i+1,flowz,id);
        }
    }
    int siz=relRoad.size();
    char topo_file[1000];
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

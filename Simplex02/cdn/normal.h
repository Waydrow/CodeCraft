vector<int> curRoad;
vector<vector<int> >relRoad;
vector<int>mustChoose;
set<int>notChoose;
map<unsigned long long,int>mp;

// char topxo_file[20000];

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


unsigned long long getHash(int *gene) {
    char hashStr[BITSIZE];
    int pos=0;
    for(int i=0; i<nodesNum; i++) {
        hashStr[pos++]=(char)(gene[i]+'a');
    }
    unsigned long long seed = 131; // 31 131 1313 13131 131313 etc..
    unsigned long long hash = 0;
    char *str=hashStr;
    while (*str) {
        hash = hash * seed + (*str++);
    }
    return hash;
}

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
        addedge(n_node,n_client,0, need);
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
vector<pair <int,int> > checkSatisfy(int nodesNum,int clientNum) {
    vector<pair<int,int> >rrr;
    memset(imNotOk,0,sizeof(imNotOk));
    int numAdd=1+nodesNum;
    //int rel=0;
    for(int i=0; i<clientNum; i++) {
        int u=i+numAdd;
        for(edge *k=V[u]; k; k=k->next) {
            if((k->t==T)&&(k->c!=0)){
                rrr.push_back(make_pair(i,k->c));
                imNotOk[u]=true;
            }
        }
    }
    return rrr;
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
void pprint(Individual &tmp){
    puts("===");
    for(int  i=0;i<nodesNum;i++){
        printf("%d ",tmp.gen[i]);
    }
    puts("");
}
void getBetter(Individual &rel) {
    int numAdd=1+nodesNum;
    for(int i=0; i<clientNum; i++) {
        int u=i+numAdd;
        int neibornode=myNodes[i];
        for(edge *k=V[u]; k; k=k->next) {
            if((k->t==T)&&(k->c!=0)){
                for(int j=0;j<serverLevel;j++){
                    if( (rel.gen[neibornode]==0&&serverInfo[j][0]>=(k->c)) ||
                        ( rel.gen[neibornode]>0&&serverInfo[j][0]-serverInfo[ rel.gen[neibornode]-1 ][0]>=(k->c) ) ){
                        rel.gen[ myNodes[i] ]=j+1;
                        break;
                    }
                    if(j==serverLevel-1){
                        int tmpflow=k->c;
                        int cha=0;
                        if(rel.gen[neibornode]){
                            cha=serverInfo[ rel.gen[neibornode]-1 ][0];
                        }
                        tmpflow-=(serverInfo[j][0]-cha);
                        rel.gen[ myNodes[i] ]=j+1;
                        int nowu=myNodes[i]+1;
                        for(edge *kk=V[nowu];kk;kk=kk->next){
                            if(tmpflow<=0)break;
                            if(kk->t>0&&kk->t<=nodesNum+1){
                                int vv=kk->t;
                                int fflow=0;
                                for(edge *gg=V[vv];gg;gg=gg->next){
                                    if(gg->t==nowu&&gg->bf>0){
                                        fflow=gg->c;
                                        break;
                                    }
                                }
                                vv--;     //二级节点编号
                                cha=0;
                                if(rel.gen[vv]){
                                    cha=serverInfo[rel.gen[vv]-1][0];
                                }
                                for(int w=0;w<serverLevel;w++){
                                    if(serverInfo[w][0]-cha >= min(tmpflow,fflow)){
                                        rel.gen[vv]=w+1;
                                        tmpflow-=(serverInfo[w][0]-cha);
                                    }
                                }
                            }
                        }
                        while(tmpflow>=0){
                            printf("%d\n",tmpflow);
                        }
                    }
                }
            }
        }
    }
}
void checkkk(){
    int numAdd=1+nodesNum;
    for(int i=0;i<clientNum;i++){
        int u=i+numAdd;
        for(edge *k=V[u]; k; k=k->next) {
            if((k->t==T)&&(k->c!=0)){
                while(1)puts("ggggg");
            }
        }
    }

}
long long mynum=0;
long long mytime=0;
int calCost(Individual gene,int mustCal,bool deleteExtraCost) {
    unsigned long long relHash=getHash(gene.gen);
    if((!mustCal)&&mp.find(relHash)!=mp.end())return mp[relHash];
    buildBasicGraph();
    int numAdd=1;

    for(int i=0; i<nodesNum; i++) {
        if(gene.gen[i]) {
            setCapacity(S,i+numAdd,serverInfo[ gene.gen[i]-1 ][0]);
        } else {
            setCapacity(S,i+numAdd,0);
        }
    }
    int costRel=maxcostflow();
    vector< pair<int,int> > poverty=checkSatisfy(nodesNum,clientNum);
    int siz=poverty.size();
    for(int i=0;i<siz;i++){
        costRel+=nodesMoney[ myNodes[ poverty[i].first ] ];
        for(int j=0;j<serverLevel;j++){
            if(serverInfo[j][0]>=poverty[i].second){
                costRel+=serverInfo[j][1];
                break;
            }
        }
    }
    for(int i=0;i<nodesNum;i++){
        if(gene.gen[i]){
            costRel+=nodesMoney[i];
            costRel+=serverInfo[ gene.gen[i]-1 ][1];
        }
    }
    return mp[relHash]=costRel;
}

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

void printRel(int *rel) {
    for(int i=0; i<nodesNum; i++) {
        if(rel[i]) {
            findRoad(i+1,serverInfo[ rel[i]-1 ][0],rel[i]-1);
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

void printGene(int *gene) {
    for(int i=0; i<nodesNum; i++) {
        cout << gene[i]<<' ';
    }
    cout << endl;
}

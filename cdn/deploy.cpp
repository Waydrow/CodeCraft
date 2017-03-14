#include "deploy.h"

#include <stdio.h>

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <queue>
#include <bitset>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <set>
#include <algorithm>
#include <map>
#include <utility>

//..xt:增广路跑费用流[SPFA]（可heap优化dij优化）
//..xt:每次需要重新建图 （可memcpy优化）
//..xt:需要保证源点编号最小，汇点编号最大

//..xt:定义最大点数边数
#define MAXN 2000
#define MAXM (52000*8)
#define INF 10000000
#define BITSIZE 1000


using namespace std;

int nodesNum; // 网络节点数量
int linkNum; // 网络链路数量
int clientNum; // 消费节点数量
int deployCost; // 服务器部署成本
vector<int> curRoad;//寻找路径dfs时记录
vector<vector<int> >relRoad;//输出答案，记录所有路径
vector<int>mustChoose;//记录一定需要选择的点的编号
map<unsigned long long,int>mp;//记录已经计算过的答案
// 最终需要输出到文件的内容
char topo_file[20000];

string relStr;

#include "population.h"

//..xt:将输入内容从topo 读入数据结构


//..xt:费用流用到全局变量以及结构体
bool vis[MAXN];
struct edge {
    edge *next,*op;
    int t,c,v,bf;
} ES[MAXM],*V[MAXN];
int N,M,S,T,EC=-1;
int demond[MAXN],sp[MAXN],Prev[MAXN];
edge *path[MAXN];
//记录链路信息
int links[50500][4];
int consumptionNodes[505][3];
int myNodes[1000];


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
        //printf("%d %d %d %d\n",st,ed,bandW,cost);
        //st+=numAdd;ed+=numAdd;
        //addedge(st,ed,cost,bandW);
        //addedge(ed,st,cost,bandW);
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
        //printf("%d %d %d\n",n_client,n_node,need);
        //n_client+=numAdd;
        //addedge(n_node,n_client,0);
        //addedge(n_client,T,0,need);
    }
}

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
    //printf("%d %d %d %d\n",a,b,v,c);
    edge e1= {V[a],0,b,c,v,c},e2= {V[b],0,a,0,-v,0};
    ES[++EC]=e1;
    V[a]=&ES[EC];
    ES[++EC]=e2;
    V[b]=&ES[EC];
    V[a]->op=V[b];
    V[b]->op=V[a];
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
    //for(int i=0;i<nodesNum;i++){
    //addedge(S,i+numAdd,0);
    //}
    for(int i=0; i<nodesNum; i++) {
        addedge(S,i+1,0,0);
    }
    for(int i=0; i<linkNum; i++) {
        //sscanf(topo[i+lineAdd],"%d%d%d%d",&st,&ed,&bandW,&cost);

        st=links[i][0];
        ed=links[i][1];
        bandW=links[i][2];
        cost=links[i][3];
        //printf("%d %d %d %d\n",st,ed,bandW,cost);
        st+=numAdd;
        ed+=numAdd;
        addedge(st,ed,cost,bandW);
        addedge(ed,st,cost,bandW);

    }
    numAdd+=nodesNum;
    int n_client,n_node,need;
    for(int i=0; i<clientNum; i++) {
        //sscanf(topo[i+linkNum],"%d%d%d",&n_client,&n_node,&need);

        n_client=consumptionNodes[i][0];
        n_node=consumptionNodes[i][1];
        need=consumptionNodes[i][2];
        //printf("%d %d %d\n",n_client,n_node,need);
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








//..xt：寻找最短路并且记录路径
//..xt:如果存在从源点到汇点大路径则返回true，否则返回false
bool SPFA() {
    int u,v;
    for(u=S; u<=T; u++) {
        sp[u]=INF;
    }
    queue<int>q;
    Prev[S]=-1;
    q.push(S);
    sp[S]=0;
    vis[S]=1;
    while(!q.empty()) {
        u=q.front();
        vis[u]=0;
        q.pop();
        for(edge *k=V[u]; k; k=k->next) {
            v=k->t;
            if(k->c>0&&sp[u]+k->v<sp[v]) {
                sp[v]=sp[u]+k->v;
                Prev[v]=u;
                path[v]=k;
                if(vis[v]==0) {
                    vis[v]=1;
                    q.push(v);
                }
            }
        }
    }
    return sp[T]!=INF;
}



//..xt:返回流最大的情况下的最小费用
int argument() {
    int i,cost=INF,flow=0;
    edge *e;
    for(i=T; Prev[i]!=-1; i=Prev[i]) {
        e=path[i];
        if(e->c<cost)cost=e->c;
    }
    //cout << cost << endl;
    for(int i=T; Prev[i]!=-1; i=Prev[i]) {
        e=path[i];
        e->c-=cost;
        e->op->c+=cost;
        //cout << e->v << endl;
        flow+=e->v*cost;
    }
    //cout << flow << endl;
    return flow;
}

//..xt:计算最大流
int maxcostflow() {
    int Flow=0;
    while(SPFA()) {
        //puts("spfa");
        //printf("%d %d\n",S,T);
        Flow+=argument();
    }
    return Flow;
}


//..xt:检测是不是所有消费节点都得到满足
int checkSatisfy(int nodesNum,int clientNum) {
    int numAdd=1+nodesNum;
    int rel=0;
    for(int i=0; i<clientNum; i++) {
        int u=i+numAdd;
        for(edge *k=V[u]; k; k=k->next) {
            if((k->t==T)&&(k->c!=0)) rel++;
        }
    }
    return rel;
}



//..xt:输出结果之前跑一此进行优化
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


//..xt:给出某个体的评估值
int calCost(bitset<BITSIZE>gene,int mustCal) {
    //printf("%d %d %d\n",nodesNum,linkNum,clientNum);
    /*for(int i=0;i<nodesNum;i++){
      if(gene[i])printf("1");
      else printf("0");
    }
    puts("");*/
    //exit(0);
    unsigned long long relHash=getHash(gene);
    if((!mustCal)&&mp.find(relHash)!=mp.end())return mp[relHash];
    buildBasicGraph(nodesNum,linkNum,clientNum);
    //puts("buildok");
    int numAdd=1;
    for(int i=0; i<nodesNum; i++) {
        //exit(0);
        if(gene[i]) {
            setCapacity(S,i+numAdd,INF);
        } else {
            setCapacity(S,i+numAdd,0);
        }
    }
    //puts("graphok");
    int costRel=maxcostflow();
    //puts("costok");
    //printf("%d\n",costRel);
    //printf("%d\n",checkSatisfy(nodesNum,clientNum));
    //exit(0);
    int addit=checkSatisfy(nodesNum,clientNum);
    /*
    if(!mustCal)
        return mp[relHash]=costRel+(gene.count())*deployCost+addit*deployCost/2;
    else */
    return mp[relHash]=costRel+(gene.count()+addit)*deployCost;
}


//你要完成的功能总入口
/* topo[MAX_EDGE_NUM] 中每一项存储输入文件中的一行
   line_num 为输入文件中的总行数
   filename 为输出文件名
*/

//..xt:递归寻找最佳解
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

//..xt：打印基因用于调试
void printGene(bitset<BITSIZE>gene) {
    for(int i=0; i<nodesNum; i++) {
        cout << gene[i];
    }
    cout << endl;
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
    //cout<<siz<<endl<<endl;
    for(int i=0; i<siz; i++) {
        int mm=relRoad[i].size();
        for(int j=0; j<mm; j++) {
            //printf("%d ",relRoad[i][j]);
            sprintf(topo_file, "%d ", relRoad[i][j]);
            relStr+=string(topo_file);
        }
        if(i<siz-1) {
            sprintf(topo_file, "\n");
            relStr+=string(topo_file);
        }
        //puts("");
    }
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
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename) {



    sscanf(topo[0], "%d%d%d", &nodesNum, &linkNum, &clientNum);
    sscanf(topo[2], "%d", &deployCost);
    readData(topo,nodesNum,linkNum,clientNum);
    getMustChoose();
    /*
    cout<< "vector:\n";
    for(int i=0;i<mustChoose.size();i++){
        printf("%d ",mustChoose[i]);
    }
    puts("");
    */
    //cout<<"readdata success"<<endl;
    Population p = Population();

    p.epoch();
    //p.show();
    bitset<BITSIZE>rel = p.everBestIndividual.bitIn;

    calCost(rel,1);
    //cout <<calCost(rel)<<endl;
    rel|=getBetter();
    calCost(rel,1);
    //cout<<"服务器个数: "<<rel.count()<<endl;
    //cout <<calCost(rel)<<endl;
    printRel(rel);
    //cout<<relStr<<endl;
    //printf("%d %d %d\n", nodesNum, linkNum, clientNum);
    //printf("%d\n", deployCost);
    //readData(topo,nodesNum,linkNum,clientNum);
    //buildBasicGraph(nodesNum,linkNum,clientNum);

    // 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
    write_result(relStr.c_str(), filename);

}

#ifndef _MinCostFlowSolution
#define _MinCostFlowSolution



#include <iostream>
#include <limits>
#include <bitset>
#include <map>
#include <utility>
#include <string.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <set>

#define BITSIZE 1300
using namespace std;
struct sst {
    int id,dis;
    sst(int a,int b) {
        id=a;
        dis=b;
    }
    sst() {};
};
struct cmp {
    bool operator()(const sst &a,const sst &b) {
        return a.dis>b.dis;
    }
};
struct eedge {
    int id,w,cost;
    eedge *next;
};


class MinCostFlowSolution {

public:

    double sumtime,times;
    void findOkPath(int);
    int findTheRoad(int,int);
    void getRestCost();
    void solutionMinCostFlow();

    void getArcFlow( double * F , int* nms = NULL ,
                  const int strt = 0 , int stp = Inf<int>() );

    inline int getArcStartNode( const int i );

    inline int getArcEndNode( const int i );

    void getArcsCost( double * Costv , const int * nms = NULL ,
                   const int strt = 0 , int stp = Inf<int>() );

    inline double getArcCost( const int i );
    template <typename T>
    class Inf {
    public:
        Inf() {}
        operator T() {
            return( std::numeric_limits<T>::max() );
        }
    };

    MinCostFlowSolution( const int nmx = 0 , const int mmx = 0 ) {
        sumtime=0;
        times=0;
        n = nmx;
        m = mmx;
        n = m = 0;
        status = -1;

        mmmin = numeric_limits<double>::epsilon()*100;
        g_is_ifirst = 1;


        nodesP = NULL;
        arcAddress = NULL;
        selectAddress = NULL;

        if( n && m )
            assignSpace();
        else
            n = m = 0;
    }
//    void getNotChoose(set<int>&notChoose);
    inline void LoadDMX( bitset<BITSIZE> &gens);
    inline unsigned long long getHash(bitset<BITSIZE>&gene);
    inline void ReadData(char *topo[], int nodesNum, int linkNum, int clientNum);
    inline pair<pair<long long,long long>, bitset<BITSIZE> > CalCost(bitset<BITSIZE> gens, int mustCal,bool deleteExtraCost);
    inline void addeedge(int a,int b,int c,int d);
    inline void findShortestPath(int num);
    inline void buidShortestMap();
    inline void getClassified();
    inline int uniform_iint(int a, int b);
    inline bitset<BITSIZE> initialBetter();
    ~MinCostFlowSolution() {
        printf("time:%.6lf\n",sumtime/times);
        deleteSelectSpace();
        deleteSpace();
    }

protected:

    template<class T>
    inline bool ETZ( T x , const T eps ) {
        if( numeric_limits<T>::is_integer )
            return( x == 0 );
        else
            return( (x <= eps ) && ( x >= -eps ) );
    }

    template<class T>
    inline bool GTZ( T x , const T eps ) {
        if( numeric_limits<T>::is_integer )
            return( x > 0 );
        else
            return( x > eps );
    }

    template<class T>
    inline bool GEZ( T x , const T eps ) {
        if( numeric_limits<T>::is_integer )
            return( x >= 0 );
        else
            return( x >= - eps );
    }

    template<class T>
    inline bool LTZ( T x , const T eps ) {
        if( numeric_limits<T>::is_integer )
            return( x < 0 );
        else
            return( x < - eps );
    }

    template<class T>
    inline bool LEZ( T x , const T eps ) {
        if( numeric_limits<T>::is_integer )
            return( x <= 0 );
        else
            return( x <= eps );
    }

    vector<int> NodesNeighbor[BITSIZE];
    vector<int> roadCost[BITSIZE];
    vector<int> roadBanw[BITSIZE];
    //vector<int,double> myVotes[BITSIZE];
    int serverLevel;//服务器类型数量
    int n;      //点数
    bool imHung[BITSIZE];
    int m;
    int restCost;
    int sumCost;
    double mmmin;
    vector<int> classf[BITSIZE][5];
    int ednum;
    eedge *aadj[1000010];
    eedge eedges[1000010];
    int dis[1000010];
    int mynum[BITSIZE];
    int nnn;
    const int inf=19999999;
    int dis1[1000050];
    int links[50500][4];
    int consumptionNodes[505][3];
    int myNodes[BITSIZE];
    int nodesNum;
    int linkNum;
    int clientNum;
    int allNeed;
    int deployCost;

    int g_is_ifirst;

    map<unsigned long long,int>mp;

    int status;

private:
    struct arcInfo;
    struct nodeInfo {
        nodeInfo *prevInT;
        nodeInfo *nextInT;
        arcInfo *enteringTArc;
        double balance;
        double potential;
        int subTreeLevel;
    };
    struct arcInfo {
        nodeInfo *wei;
        nodeInfo *tou;
        double flow;
        double cost;

        char flagInTree;
        double upper;
    };
    struct arcSelect {
        arcInfo *arc;

        double val;           // 调整后减少费用的绝对值
    };
    nodeInfo *nodesP;       //指向点的链表，指向root人工点
    nodeInfo *humanNode;
    nodeInfo *illegalNodeAddress;         // 第一个点的非法地址
    arcInfo *arcAddress;               // 链表头，指向最后一条边
    arcInfo *humanArc;          // 链表头，指向人工边的最后一条
    arcInfo *illegalArcAddress;           // 第一个边的非法地址
    arcInfo *illegalDummyArcAddress;
    arcInfo *beginArc;
    double iterator;
    arcSelect *selectAddress;       // 限制调整的边
    int amountOfGroup;
    int tempSelectSize;
    int positionOfGroup;
    int amountSelectList;
    int popularSize;
    double *balAddress;


    void assignSpace();
    void deleteSpace(); //删除边链表，点链表，候选边链表
    void assignSelectSpace();
    void deleteSelectSpace();
    void initialAlg();
    void simplex();

    template<class N, class A>
    void getNewTree( A *h , A *k , N *h1 , N *h2 , N *k1 , N *k2 );

    template<class N>
    N* solveTree( N *root, int delta );

    template<class N>
    void linkTree( N *root , N *lastNode , N *previousNode );

    arcInfo* getEnteringArc( void );
    inline void generateSelect( void );
    inline void orderSelect( int min , int max );

    template<class N, class RCT>
    inline void changeP( N *r , RCT delta );

    template<class N, class A>
    inline N* findF( N *n, A *a );

    template<class A>
    inline double downCost( A *a );
};


inline int MinCostFlowSolution::getArcStartNode( const int i ) {
    return( int( ( (arcAddress + i)->wei - nodesP + 1 ) ) );
}

inline int MinCostFlowSolution::getArcEndNode( const int i ) {
    return( int( ( (arcAddress + i)->tou - nodesP + 1 )  ) );
}

inline double MinCostFlowSolution::getArcCost( const int i ) {
    return( (arcAddress + i)->cost );
}

template <class A>
inline double MinCostFlowSolution::downCost( A *a ) {
    double redc = (a->wei)->potential - (a->tou)->potential;
    redc = redc + a->cost;
    return( redc );
}

inline bitset<BITSIZE> MinCostFlowSolution::initialBetter(){
    bitset<BITSIZE>gene;
    gene.reset();
    for(int i=0;i<clientNum;i++){
        int nextPos=uniform_iint(0,2);
        for(int j=nextPos;;){
            if(classf[i][j].size()>0){
                int another=uniform_iint(0,((int)classf[i][j].size())-1);
                int relPos=classf[i][j][another];
                gene[relPos]=1;
                break;
            }
            j++;
            if(j>=3)j-=3;
            if(j==nextPos)break;
        }
    }
    return gene;
}

inline void MinCostFlowSolution::addeedge(int a,int b,int c,int d=0) {
    eedge *aa;
    aa=&eedges[ednum];
    ednum++;
    aa->id=b;
    aa->w=c;
    aa->next=aadj[a];
    aadj[a]=aa;
}

inline void MinCostFlowSolution::findShortestPath(int num) {
    num--;
    for(int i=0; i<nnn; i++) {
        dis[i]=inf;
    }
    priority_queue<sst,vector<sst>,cmp>q;
    sst tmp;
    dis[num]=0;
    q.push(sst(num,0));
    while(!q.empty()) {
        tmp=q.top();
        q.pop();
        for(eedge *p=aadj[tmp.id]; p; p=p->next) {
            if(dis[p->id]>p->w+dis[tmp.id]) {
                dis[p->id]=p->w+dis[tmp.id];
                q.push(sst(p->id,dis[p->id]));
            }
        }
    }
    for(int i=0; i<nodesNum; i++) {
        deployCost = serverInfo[serverLevel-1][1] + nodesMoney[ myNodes[num-nodesNum] ];
        double kkk=dis[i]*consumptionNodes[num-nodesNum][2];
        kkk += nodesMoney[i]  ;
        if(kkk>=deployCost)continue;
        double gailv=(deployCost-kkk)/deployCost;
        if(gailv>=0.7)classf[num-nodesNum][0].push_back(i);
        else if(gailv>=0.3)classf[num-nodesNum][1].push_back(i);
        else classf[num-nodesNum][2].push_back(i);
    }
}

inline void MinCostFlowSolution::buidShortestMap() {
    ednum=0;
    memset(aadj,0,sizeof(aadj));
    for(int i=0; i<linkNum; i++) {
        addeedge(links[i][0],links[i][1],links[i][3]);
        addeedge(links[i][1],links[i][0],links[i][3]);
    }
    for(int i=0; i<clientNum; i++) {
        addeedge(consumptionNodes[i][0]+nodesNum,consumptionNodes[i][1],0);
        addeedge(consumptionNodes[i][1],consumptionNodes[i][0]+nodesNum,0);
    }
    nnn=nodesNum+clientNum;

}

inline void MinCostFlowSolution::getClassified() {
    buidShortestMap();
    for(int i=nodesNum+1; i<=nodesNum+clientNum; i++) {
        findShortestPath(i);
    }
}

inline int MinCostFlowSolution::uniform_iint(int a, int b) {
    if (g_is_ifirst) {
        g_is_ifirst = 0;
        srand((unsigned int)time(NULL));
    }

    return (int)((double)rand() / ((RAND_MAX + 1.0) / (b - a + 1.0)) + a);
}



//..xt:get hash of each gen
inline unsigned long long MinCostFlowSolution::getHash(bitset<BITSIZE>&gene) {
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

inline void MinCostFlowSolution::ReadData(char *topo[], int nodesNums, int linkNums, int clientNums) {
  serverLevel=0;
  nodesNum = nodesNums;
  linkNum = linkNums;
  clientNum = clientNums;
  allNeed=0;
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
          allNeed+=tmpc;
          consumptionNodesI++;
      }
  }

    getClassified();
}


inline pair<pair<long long,long long>,bitset<BITSIZE> > MinCostFlowSolution::CalCost(bitset<BITSIZE> gens, int mustCal,bool deleteExtraCost) {

    unsigned long long relHash=getHash(gens);
    if((!mustCal)&&mp.find(relHash)!=mp.end())return make_pair(make_pair(mp[relHash],mp[relHash]),gens);
    LoadDMX(gens);
    int a=clock();
    solutionMinCostFlow();
    double b=clock()-a;
    sumtime+=b/CLOCKS_PER_SEC;
    times++;
    double * x = new double[ m ];
    double * c = new double[ m ];
    getArcsCost(c);
    getArcFlow(x);
    int myflow[BITSIZE];
    memset(myflow,0,sizeof(myflow));
    double sum = 0;

    int myOutSum[BITSIZE];
    int myInSum[BITSIZE];
    memset(myOutSum,0,sizeof(myOutSum));
    memset(myInSum,0,sizeof(myInSum));
    for( int i = 0 ; i < m ; i++ ) {
        int sst=getArcStartNode(i);
        int eed=getArcEndNode(i);
        if(sst>1&&sst<nodesNum+2){
            myflow[sst-2]+=x[i];
            if(eed>1&&eed<nodesNum+2){
                myflow[eed-2]-=x[i];
            }
        }
        if(eed>nodesNum+1) {
            myInSum[eed-nodesNum-2] += x[i];
        }
        if(sst>1&&sst<nodesNum+2&&gens[sst-2]) {
            myOutSum[sst-2] += x[i];
        }
        if (c[i] != 4e7) {
            sum += c[i] * x[i];
        }
        if(getArcStartNode(i)==1 && getArcEndNode(i)>nodesNum+1 && x[i]>0) {
            if(gens[ myNodes[eed-nodesNum-2] ]){
                while(true){
                    puts("adderror1");
                }
            }
            for(int w=0;w<serverLevel;w++){
                if(serverInfo[w][0]>=x[i]){
                    sum+=serverInfo[w][1];
                    sum+=nodesMoney[ myNodes[eed-nodesNum-2] ];
                    break;
                }
                if(w==serverLevel-1){
                    while(true){
                        puts("adderror2");
                    }
                }
            }
        }

    }
    //getRestCost();
    //sum-=restCost
    for(int i=0;i<nodesNum;i++){
        while(myflow[i]<0 || myflow[i]>serverInfo[serverLevel-1][0]){
            puts("aaaflowadderror");
        }
        if(gens[i]){
            sum+=nodesMoney[i];
            for(int w=0;w<serverLevel;w++){
                if(serverInfo[w][0]>=myflow[i]){
                    sum+=serverInfo[w][1];
                    break;
                }
                while(w==serverLevel-1){
                    puts("adderror3");
                }
            }
        }
    }
    long long tmpsum=sum;
    int mmin=4e7,outPos=0;
    for(int i=0; i<nodesNum; i++) {
        if(gens[i]) {
            if(myOutSum[i]==0) {
                gens[i]=0;
            }
            if(mmin>myOutSum[i]) {
                mmin=myOutSum[i];
                outPos=i;
            }
        }
    }
    if(mmin<500){
        for(int w=0;w<serverLevel;w++){
            if(serverInfo[w][0]>=myOutSum[outPos]){
                tmpsum-=serverInfo[w][1];;
                tmpsum-=nodesMoney[outPos];
                break;
            }
        }
        gens[outPos]=0;
    }
    mmin=-1;
    for(int i=0; i<clientNum; i++) {
        if(mmin<consumptionNodes[i][2]-myInSum[i]) {
            mmin=consumptionNodes[i][2]-myInSum[i];
            outPos=i;
        }
    }

    if(mmin>0){
        int nextPos=uniform_iint(0,2);
        for(int i=nextPos;;){
            if(classf[outPos][i].size()>0){
                int another=uniform_iint(0,((int)classf[outPos][i].size())-1);
                int relPos=classf[outPos][i][another];
                tmpsum+=serverInfo[serverLevel-1][1];
                tmpsum+=nodesMoney[relPos];
                gens[relPos]=1;
                break;
            }
            i++;
            if(i>=3)i-=3;
            if(i==nextPos)break;
        }
    }
    delete[] x;
    delete[] c;
    mp[relHash] = sum;
    return make_pair(make_pair((long long)sum,tmpsum),gens);
}

inline void MinCostFlowSolution::LoadDMX( bitset<BITSIZE> &gens) {

    int serverNum = gens.count();
    int tn = nodesNum + clientNum + 1;

    int tm = linkNum*2 + serverNum + clientNum*2;

    deleteSelectSpace();
    if( n && m )  { //如果原先已经建图，删除内存
        deleteSpace();
        n = m = 0;
    }
    //申请新的内存
    n = tn;
    m = tm;
    assignSpace();
    //初始化点数和边数
    illegalNodeAddress = nodesP + n;
    humanNode = nodesP + n;

    //初始化需求量
    nodeInfo *node = nodesP;
    node->balance= -allNeed;
    node++;
    for(int i=1;i<=nodesNum;i++){
        node->balance=0;
        node++;
    }
    for(int j=nodesNum+2;j<=nodesNum+clientNum+1;j++){
        node->balance = consumptionNodes[j-(2+nodesNum)][2];
        node++;
    }

    illegalArcAddress = arcAddress + m;
    humanArc = arcAddress + m;
    illegalDummyArcAddress = humanArc + n;
    arcInfo *arc = arcAddress ;
    /*建图过程*/
    for (int j = 0; j<nodesNum; j++) {
        if (gens[j]) {
            arc->wei = nodesP+1-1;
            arc->tou = nodesP+j+2-1;
            arc->upper = serverInfo[serverLevel-1][0];
            arc->cost = 0;
            arc++;
        }
    }
    for (int j = 0; j < linkNum; j++) {
        arc->wei = nodesP+links[j][0]+2-1;
        arc->tou = nodesP+links[j][1]+2-1;
        arc->upper = links[j][2];
        arc->cost = links[j][3];
        arc++;
        arc->wei = nodesP+links[j][1]+2-1;
        arc->tou = nodesP+links[j][0]+2-1;
        arc->upper = links[j][2];
        arc->cost = links[j][3];
        arc++;
    }
    for (int j = 0; j < clientNum; j++) {
        arc->wei = nodesP+consumptionNodes[j][1]+2-1;
        arc->tou = nodesP+consumptionNodes[j][0] + 2 + nodesNum-1;
        arc->upper = consumptionNodes[j][2];
        arc->cost = 0;
        arc++;
        arc->wei = nodesP+1-1;
        arc->tou = nodesP+consumptionNodes[j][0]+2+nodesNum-1;
        arc->upper = consumptionNodes[j][2];
        arc->cost = 40000000;
        arc++;
    }

    assignSelectSpace();
    status = 0;
}

template<class T>
inline T ABS( const T x ) {
    return( x >= T( 0 ) ? x : - x );
}

template<class T>
inline void Swap( T &v1 , T &v2 ) {
    T temp = v1;
    v1 = v2;
    v2 = temp;
}


void MinCostFlowSolution::findOkPath(int num){
  for(int i=0; i<nnn; i++) {
      dis[i]=inf;
  }
  priority_queue<sst,vector<sst>,cmp>q;
  sst tmp;
  dis[num]=0;
  q.push(sst(num,0));
  while(!q.empty()) {
      tmp=q.top();
      q.pop();
      for(eedge *p=aadj[tmp.id]; p; p=p->next) {
          if(dis[p->id]>p->w+dis[tmp.id]) {
              dis[p->id]=p->w+dis[tmp.id];
              q.push(sst(p->id,dis[p->id]));
          }
      }
  }
  for(int i=1; i<=nodesNum; i++) {
      double kkk=dis[i]*consumptionNodes[num-nodesNum-1][2];
      if(kkk>=deployCost)continue;
      double gailv=(deployCost-kkk)/deployCost;
      if(gailv<0.5){
          mynum[i-1]++;
      }
  }
}
/*
void MinCostFlowSolution::getNotChoose(set<int>&notChoose){
    buidShortestMap();
    memset(mynum,0,sizeof(mynum));
    for(int i=nodesNum+1; i<=nodesNum+clientNum; i++) {
        findOkPath(i);
    }
    int ave=0;
    for(int i=0;i<nodesNum;i++){
        ave+=mynum[i];
    }
    ave/=nodesNum;
    ave/=1.15;
    for(int i=0;i<nodesNum;i++){
        //printf("mynum[%d]=%d\n",i,mynum[i]);
        if(mynum[i]<ave){
            notChoose.insert(i);
        }
    }
    cout << nodesNum << endl;
    cout << notChoose.size() <<endl;
    //exit(0);
}
*/
int MinCostFlowSolution::findTheRoad(int pos,int flow){
    if(pos>nodesNum){
        if(imHung[pos]){
            restCost+=flow*sumCost;
        }
        return flow;
    }
    int totalFlow=0;
    for(eedge *i=aadj[pos];i;i=i->next){
        if(i->w > 0){
            sumCost+=i->cost;
            int tmpFlow=findTheRoad(i->id,min(flow,i->w));
            sumCost-=i->cost;
            totalFlow+=tmpFlow;
            i->w -= tmpFlow;
        }
    }
    return totalFlow;
}

void MinCostFlowSolution::getRestCost(){
    restCost=0;
    sumCost=0;
    ednum=0;
    memset(imHung,0,sizeof(imHung));
    memset(aadj,0,sizeof(aadj));
    for( arcInfo *arc = arcAddress ; arc != illegalArcAddress ; arc++ ){
        int sst=arc->tou - nodesP;
        int eed=arc->wei - nodesP;
        int x=arc->flow;
        int cost=arc->cost;
        if(sst==0&&eed>nodesNum&&x>0){
            imHung[eed]=true;
            continue;
        }
        if(sst>nodesNum+clientNum||eed>nodesNum+clientNum){
            continue;
        }
        addeedge(sst,eed,x,cost);
    }
    findTheRoad(0,100000000);
}

void MinCostFlowSolution::solutionMinCostFlow( void ) {
    if( status == 0 ){
        initialAlg();
    }
    simplex();
}

void MinCostFlowSolution::getArcFlow( double * F , int* nms , const int strt , int stp ) {
    if( stp > m )
        stp = m;

    if( nms ) {
        for( int i = strt ; i < stp ; i++ ) {
            double tXi = ( arcAddress + i )->flow;
            if( GTZ( tXi , mmmin ) ) {
                *(F++) = tXi;
                *(nms++) = i;
            }
        }
        *nms = Inf<int>();
    }else
        for( int i = strt; i < stp; i++ )
            *(F++) = ( arcAddress + i )->flow;

}

void MinCostFlowSolution::getArcsCost( double * Costv , const int * nms ,
                           const int strt , int stp ) {
    if( stp > m )
        stp = m;

    if( nms ) {
        while( *nms < strt )
            nms++;

        for( int h ; ( h = *(nms++) ) < stp ; )
            *(Costv++) = (arcAddress + h)->cost;
    } else
        for( arcInfo* arc = arcAddress + strt ; arc < (arcAddress + stp) ; arc++ )
            *(Costv++) = arc->cost;

}

void MinCostFlowSolution::assignSpace( void ) {
    nodesP = new nodeInfo[ n + 1 ];
    arcAddress = new arcInfo[ m + n ];
    humanArc = arcAddress + m;
}

void MinCostFlowSolution::deleteSpace() {
    delete[] nodesP;
    delete[] arcAddress;
    nodesP = NULL;
    arcAddress = NULL;
    deleteSelectSpace( );
}

void MinCostFlowSolution::assignSelectSpace() {
    if( m < 10000 ) {
        amountSelectList = 30;
        popularSize = 5;
    } else if( m > 100000 ) {
        amountSelectList = 200;
        popularSize = 20 ;
    } else {
        amountSelectList = 50;
        popularSize = 10;
    }
    selectAddress = new arcSelect[ popularSize + amountSelectList + 1 ];
}

void MinCostFlowSolution::deleteSelectSpace( void ) {
    delete[] selectAddress;
    selectAddress = NULL;
}

void MinCostFlowSolution::initialAlg( void ) {
    arcInfo *arc;
    nodeInfo *node;
    for( arc = arcAddress ; arc != illegalArcAddress ; arc++ ) {
        arc->flow = 0;
        arc->flagInTree = 1;
    }

    for( arc = humanArc ; arc != illegalDummyArcAddress ; arc++ ) {
        node = nodesP + ( arc - humanArc );
        if( node->balance > 0 ) {
            arc->wei = humanNode;
            arc->tou = node;
            arc->flow = node->balance;
        } else {
            arc->wei = node;
            arc->tou = humanNode;
            arc->flow = -node->balance;
        }

        arc->cost = (1e10);
        arc->flagInTree = 0;
        arc->upper = arc->flow;
    }

    humanNode->balance = 0;
    humanNode->prevInT = NULL;
    humanNode->nextInT = nodesP;
    humanNode->enteringTArc = NULL;
    humanNode->potential = (1e10);
    humanNode->subTreeLevel = 0;
    for( node = nodesP ; node != illegalNodeAddress ; node++) {
        node->prevInT = node - 1;
        node->nextInT = node + 1;
        node->enteringTArc = humanArc + (node - nodesP);
        if( node->balance > 0 )
            node->potential = 2 * (1e10);
        else
            node->potential = 0;

        node->subTreeLevel = 1;
    }
    nodesP->prevInT = humanNode;
    ( nodesP + n - 1 )->nextInT = NULL;
}

void MinCostFlowSolution::simplex( void ) {

    status = 0;
    beginArc = arcAddress;

    iterator = 0;
    arcInfo *enteringArc;
    arcInfo *leavingArc;
    generateSelect();
    while( status == 0 ) {
        iterator++;
        enteringArc = getEnteringArc();

        if( enteringArc ) {
            arcInfo *arc;
            nodeInfo *k1;
            nodeInfo *k2;
            double t;
            double theta;
            if( enteringArc->flagInTree == 2 ) {
                k1 = enteringArc->tou;
                k2 = enteringArc->wei;
                theta = enteringArc->flow;
            } else {
                k1 = enteringArc->wei;
                k2 = enteringArc->tou;
                theta = enteringArc->upper - enteringArc->flow;
            }
            nodeInfo *memK1 = k1;
            nodeInfo *memK2 = k2;
            leavingArc = NULL;
            bool leavingReducesFlow = GTZ( downCost( enteringArc ) , mmmin );
            bool leave;
            while( k1 != k2 ) {
                if( k1->subTreeLevel > k2->subTreeLevel ) {
                    arc = k1->enteringTArc;
                    if( arc->wei != k1 ) {
                        t = arc->upper - arc->flow;
                        leave = false;
                    } else {
                        t = arc->flow;
                        leave = true;
                    }
                    if( t < theta ) {
                        theta = t;
                        leavingArc = arc;
                        leavingReducesFlow = leave;
                    }
                    k1 = findF( k1 , arc );
                } else {
                    arc = k2->enteringTArc;
                    if( arc->wei == k2 ) {
                        t = arc->upper - arc->flow;
                        leave = false;
                    } else {
                        t = arc->flow;
                        leave = true;
                    }
                    if( t <= theta ) {
                        theta = t;
                        leavingArc = arc;
                        leavingReducesFlow = leave;
                    }
                    k2 = findF(k2, arc);
                }
            }
            if( leavingArc == NULL )
                leavingArc = enteringArc;

            k1 = memK1;
            k2 = memK2;

            if( ! ETZ(theta , mmmin ) ) {
                if( enteringArc->wei == k1 )
                    enteringArc->flow = enteringArc->flow + theta;
                else
                    enteringArc->flow = enteringArc->flow - theta;

                while( k1 != k2 ) {
                    if( k1->subTreeLevel > k2->subTreeLevel ) {
                        arc = k1->enteringTArc;
                        if( arc->wei != k1 )
                            arc->flow = arc->flow + theta;
                        else
                            arc->flow = arc->flow - theta;

                        k1 = findF(k1, k1->enteringTArc);
                    } else {
                        arc = k2->enteringTArc;
                        if( arc->wei == k2 )
                            arc->flow = arc->flow + theta;
                        else
                            arc->flow = arc->flow - theta;

                        k2 = findF(k2, k2->enteringTArc);
                    }
                }
            }

            if( enteringArc != leavingArc ) {
                bool leavingBringFlowInT2 = ( leavingReducesFlow ==
                                              ( ( leavingArc->wei )->subTreeLevel > ( leavingArc->tou )->subTreeLevel ) );
                if( leavingBringFlowInT2 == ( memK1 == enteringArc->wei ) ) {
                    k2 = enteringArc->wei;
                    k1 = enteringArc->tou;
                } else {
                    k2 = enteringArc->tou;
                    k1 = enteringArc->wei;
                }
            }
            if( leavingReducesFlow )
                leavingArc->flagInTree = 1;
            else
                leavingArc->flagInTree = 2;

            if( leavingArc != enteringArc ) {
                enteringArc->flagInTree = 0;
                nodeInfo *h1;
                nodeInfo *h2;
                if( ( leavingArc->wei )->subTreeLevel < ( leavingArc->tou )->subTreeLevel ) {
                    h1 = leavingArc->wei;
                    h2 = leavingArc->tou;
                } else {
                    h1 = leavingArc->tou;
                    h2 = leavingArc->wei;
                }

                getNewTree(leavingArc, enteringArc, h1, h2, k1, k2);
                k2 = enteringArc->tou;
                double delta = downCost(enteringArc);
                if( ( enteringArc->wei )->subTreeLevel > ( enteringArc->tou )->subTreeLevel ) {
                    delta = -delta;
                    k2 = enteringArc->wei;
                }

                changeP( k2 , delta );
            }
        } else {
            status = 1;
        }
    }

}

template<class N, class A>
void MinCostFlowSolution::getNewTree( A *h , A *k , N *h1 , N *h2 , N *k1 , N *k2 ) {
    int delta = (k1->subTreeLevel) + 1 - (k2->subTreeLevel);
    N *root = k2;
    N *dad;
    N *previousNode = k1;
    N *lastNode;
    A *arc1 = k;
    A *arc2;
    bool fine = false;
    while( fine == false ) {
        if( root == h2 )
            fine = true;

        dad = findF( root , root->enteringTArc );
        lastNode = solveTree( root , delta );
        linkTree( root , lastNode , previousNode );
        previousNode = lastNode;
        delta = delta + 2;
        arc2 = root->enteringTArc;
        root->enteringTArc = arc1;
        arc1 = arc2;
        root = dad;
    }
}


template<class N>
N* MinCostFlowSolution::solveTree( N *root , int delta ) {
    int level = root->subTreeLevel;
    N *node = root;
    while ( ( node->nextInT ) && ( ( node->nextInT )->subTreeLevel > level ) ) {
        node = node->nextInT;
        node->subTreeLevel = node->subTreeLevel + delta;
    }

    root->subTreeLevel = root->subTreeLevel + delta;
    if( root->prevInT )
        ( root->prevInT )->nextInT = node->nextInT;
    if( node->nextInT )
        ( node->nextInT )->prevInT = root->prevInT;

    return( node );
}

template<class N>
void MinCostFlowSolution::linkTree( N *root , N *lastNode , N *previousNode ) {

    N *nextNode = previousNode->nextInT;
    root->prevInT = previousNode;
    previousNode->nextInT = root;
    lastNode->nextInT = nextNode;
    if( nextNode )
        nextNode->prevInT = lastNode;
}
MinCostFlowSolution::arcInfo* MinCostFlowSolution::getEnteringArc() {
    int next = 0;
    int i;
    int minimeValue;
    if( popularSize < tempSelectSize )
        minimeValue = popularSize;
    else
        minimeValue = tempSelectSize;
    for( i = 2 ; i <= minimeValue ; i++ ) {
        arcInfo *arc = selectAddress[i].arc;
        double red_cost = downCost( arc );

        if( ( LTZ( red_cost , mmmin ) && ( arc->flagInTree == 1 ) ) ||
                ( GTZ( red_cost , mmmin ) && ( arc->flagInTree == 2 ) ) ) {
            next++;
            selectAddress[ next ].arc = arc;
            selectAddress[ next ].val = ABS( red_cost );
        }
    }

    tempSelectSize = next;
    int oldGroupPos = positionOfGroup;
    do {
        arcInfo *arc;
        for( arc = arcAddress + positionOfGroup ; arc < illegalArcAddress ; arc += amountOfGroup ) {
            if( arc->flagInTree == 1 ) {
                double red_cost = downCost( arc );
                if( LTZ( red_cost , mmmin ) ) {
                    tempSelectSize++;
                    selectAddress[ tempSelectSize ].arc = arc;
                    selectAddress[ tempSelectSize ].val = ABS( red_cost );
                }
            } else if( arc->flagInTree == 2 ) {
                double red_cost = downCost( arc );
                if( GTZ( red_cost , mmmin ) ) {
                    tempSelectSize++;
                    selectAddress[ tempSelectSize ].arc = arc;
                    selectAddress[ tempSelectSize ].val = ABS( red_cost );
                }
            }
        }

        positionOfGroup++;
        if( positionOfGroup == amountOfGroup )
            positionOfGroup = 0;

    } while( ( tempSelectSize < popularSize ) && ( positionOfGroup != oldGroupPos ) );

    if( tempSelectSize ) {
        orderSelect( 1 , tempSelectSize );
        return( selectAddress[ 1 ].arc );
    } else
        return( NULL );
}
inline void MinCostFlowSolution::generateSelect( void ) {
    amountOfGroup = ( ( m - 1 ) / amountSelectList ) + 1;
    positionOfGroup = 0;
    tempSelectSize = 0;
}

inline void MinCostFlowSolution::orderSelect( int min , int max ) {
    int left = min;
    int right = max;
    double cut = selectAddress[ ( left + right ) / 2 ].val;
    do {
        while( selectAddress[ left ].val > cut)
            left++;
        while( cut > selectAddress[ right ].val)
            right--;

        if( left < right )
            Swap( selectAddress[ left ] , selectAddress[ right ] );

        if(left <= right) {
            left++;
            right--;
        }
    } while( left <= right );

    if( min < right )
        orderSelect( min , right );
    if( ( left < max ) && ( left <= popularSize ) )
        orderSelect( left , max );
}

template<class N, class RCT>
inline void MinCostFlowSolution::changeP( N *r , RCT delta ) {
    int level = r->subTreeLevel;
    N *n = r;

    do {
        n->potential = n->potential + delta;
        n = n->nextInT;
    } while ( ( n ) && ( n->subTreeLevel > level ) );
}
template<class N, class A>
inline N* MinCostFlowSolution::findF( N *n , A *a ) {
    if( a == NULL )
        return NULL;

    if( a->wei == n )
        return( a->tou );
    else
        return( a->wei );
}

#endif

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

#define BITSIZE 1000
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
    void getNotChoose(set<int>&notChoose);
    inline void LoadDMX( bitset<BITSIZE> &gens);
    inline unsigned long long getHash(bitset<BITSIZE>&gene);
    inline void ReadData(char *topo[], int nodesNum, int linkNum, int clientNum,int deployCost);
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

    int n;      //点数
    bool imHung[1000];
    int m;
    int restCost;
    int sumCost;
    double mmmin;
    vector<int> classf[1005][5];
    int ednum;
    eedge *aadj[1000010];
    eedge eedges[1000010];
    int dis[1000010];
    int mynum[2000];
    int nnn;
    const int inf=19999999;
    int dis1[1000050];
    int links[50500][4];
    int consumptionNodes[505][3];
    int myNodes[1000];
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
        if(true)classf[num-nodesNum-1][0].push_back(i-1);
        else if(gailv>=0.1)classf[num-nodesNum-1][1].push_back(i-1);
        else classf[num-nodesNum-1][2].push_back(i-1);
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

inline void MinCostFlowSolution::ReadData(char *topo[], int nodesNums, int linkNums, int clientNums,int deployCosts) {
    nodesNum = nodesNums;
    linkNum = linkNums;
    clientNum = clientNums;
    deployCost=deployCosts;
    int lineAdd = 4;
    int startNode, endNode, upper, cost;
    for (int i = 0; i < linkNum; i++) {
        sscanf(topo[i+lineAdd], "%d%d%d%d", &startNode, &endNode, &upper, &cost);
        links[i][0] = startNode;
        links[i][1] = endNode;
        links[i][2] = upper;
        links[i][3] = cost;
    }
    lineAdd = linkNum + 5;
    allNeed = 0;
    int n_snode, n_enode, need;
    for (int i = 0; i < clientNum; i++) {
        sscanf(topo[i+lineAdd], "%d%d%d", &n_snode, &n_enode, &need);
        myNodes[n_snode] = n_enode;
        consumptionNodes[i][0] = n_snode;
        consumptionNodes[i][1] = n_enode;
        consumptionNodes[i][2] = need;
        allNeed += need;
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

    double sum = 0;
    int povertyNum=0;
    int myOutSum[BITSIZE];
    int myInSum[BITSIZE];
    memset(myOutSum,0,sizeof(myOutSum));
    memset(myInSum,0,sizeof(myInSum));
    for( int i = 0 ; i < m ; i++ ) {
        int sst=getArcStartNode(i);
        int eed=getArcEndNode(i);
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
            povertyNum++;
        }

    }
    //getRestCost();
    //sum-=restCost
    long long tmpsum=sum;
    sum+=(gens.count()+povertyNum)*deployCost;
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
    if(mmin<500)
    gens[outPos]=0;
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
                gens[relPos]=1;
                break;
            }
            i++;
            if(i>=3)i-=3;
            if(i==nextPos)break;
        }
    }
    tmpsum += (gens.count()+max(povertyNum,povertyNum))*deployCost;
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
            arc->upper = allNeed;
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

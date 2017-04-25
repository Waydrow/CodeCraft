#ifndef _MinCostFlowSolution
#define _MinCostFlowSolution

#include <iostream>
#include <limits>
#include <map>
#include <utility>
#include <string.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <set>
// #include "Individual.h"

#define BITSIZE 3700
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


class MinCostFlowSolution {

public:

    double sumtime,times;
//    void findOkPath(int);
    int findTheRoad(int,int);
    void getRestCost();
    void solutionMinCostFlow();

    inline int getArcStartNode( const int i );

    inline int getArcEndNode( const int i );

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
    inline void LoadDMX( Individual &gens);
    inline unsigned long long getHash(int *gene);
    inline void ReadData(char *topo[], int nodesNum, int linkNum, int clientNum);
    inline void CalCost(Individual &gens, int mustCal,bool deleteExtraCost);
    inline vector<int>getLastSolution(Individual &sol);
    inline void addeedge(int a,int b,int c,int d);
    //inline void findShortestPath(int num);
    inline void buidShortestMap();
    //inline void getClassified();
    inline int uniform_iint(int a, int b);
    ~MinCostFlowSolution() {
        //printf("time:%.6lf\n",sumtime/times);
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
    int mynum[2000];
    int nnn;
    const int inf=19999999;
    int dis1[1000050];
    int links[50500][4];
    int consumptionNodes[505][3];
    int myNodes[BITSIZE];
    int serverInfo[15][2];
    int nodesMoney[BITSIZE];
    int nodesNum;
    int linkNum;
    int clientNum;
    int allNeed;

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

template <class A>
inline double MinCostFlowSolution::downCost( A *a ) {
    double redc = (a->wei)->potential - (a->tou)->potential;
    redc = redc + a->cost;
    return( redc );
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
/*
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
*/

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

/*
inline void MinCostFlowSolution::getClassified() {
    buidShortestMap();
    for(int i=nodesNum+1; i<=nodesNum+clientNum; i++) {
        findShortestPath(i);
    }
}
*/

inline int MinCostFlowSolution::uniform_iint(int a, int b) {
    if (g_is_ifirst) {
        g_is_ifirst = 0;
        srand((unsigned int)time(NULL));
    }

    return (int)((double)rand() / ((RAND_MAX + 1.0) / (b - a + 1.0)) + a);
}



//..xt:get hash of each gen
inline unsigned long long MinCostFlowSolution::getHash(int *gene) {
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

    //getClassified();
}

inline void MinCostFlowSolution::CalCost(Individual &gene, int mustCal,bool deleteExtraCost) {

    unsigned long long relHash=getHash(gene.gen);
    if((!mustCal)&&mp.find(relHash)!=mp.end()){
        gene.cost=mp[relHash];
        return;
    }
    LoadDMX(gene);
    int a=clock();
    solutionMinCostFlow();
    double b=clock()-a;
    sumtime+=b/CLOCKS_PER_SEC;
    times++;
    double sum = 0;

    for( int i = 0 ; i < m ; i++ ) {
        if ((arcAddress+i)->cost != 4e7) {
            sum += (arcAddress+i)->cost * (arcAddress+i)->flow;
        }

        //如果某消费节点的需求没满足
        if(getArcStartNode(i)==1 && getArcEndNode(i)>nodesNum+1 && (arcAddress+i)->flow > 0) {
            // bool ok=false;
            //寻找最少的补足需求的服务器等级
            for(int j=0;j<serverLevel;j++){
                if(serverInfo[j][0]>=(arcAddress+i)->flow){
                    // ok=true;
                    sum+=serverInfo[j][1];
                    break;
                }
            }
            //加上场地费用
            sum+=nodesMoney[ myNodes[ getArcEndNode(i)-nodesNum-2 ] ];
        }
    }
    //加上已经布置服务器节点的服务器费用和场地费用
    for(int i=0;i<nodesNum;i++){
        if(gene.gen[i]){
            sum+=nodesMoney[i];
            sum+=serverInfo[ gene.gen[i]-1 ][1];
        }
    }
    mp[relHash] = sum;
    gene.cost=sum;
}

inline void MinCostFlowSolution::LoadDMX( Individual &gens) {

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
        if (gens.gen[j]) {
            arc->wei = nodesP+1-1;
            arc->tou = nodesP+j+2-1;
            arc->upper = serverInfo[ gens.gen[j]-1 ][0];
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

inline vector<int> MinCostFlowSolution::getLastSolution(Individual &sol){
    vector<int>rel;
    sol.updateRealPlan();
    for(int i=0;i<nodesNum;i++){
        rel.push_back(sol.gen[i]);
    }
    return rel;
}

#endif

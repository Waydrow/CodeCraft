#include "MinCostFlowSolution.h"

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

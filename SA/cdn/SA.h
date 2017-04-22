#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <stdio.h>
#include <time.h>
#include <math.h>

// 20
#define T     3000
#define EPS   1e-8
#define DELTA 0.98
#define LIMIT 1000
#define OLOOP 1000
#define ILOOP 100
#define Iterations 800

#define P_GONE 0.8
#define P_SWAP 0.8
#define P_FLIP 0.7

int may_server_num;

using namespace std;

int iterations = 0;
bitset<BITSIZE> GetNext(bitset<BITSIZE> ans) {

	if(int(ans.count()) > may_server_num) {
		double pGone = get_random_real(0, 1);
		if (pGone < P_GONE) {
			int postion = get_random_int(0, nodesNum - 1);
			if(ans[postion]){ // 为 1 则反转
				ans.flip(postion);
			}
		}
	} else {
		double pFlip = get_random_real(0, 1);
		if(pFlip < P_FLIP) {
			int postion = get_random_int(0, nodesNum - 1);
			ans.flip(postion);
		}
	}
	double pSwap = get_random_real(0, 1);
	if(pSwap < P_SWAP) {
		int x = get_random_int(0, nodesNum - 1);
		int y = get_random_int(0, nodesNum - 1);
		int temp = ans[x];
		ans[x] = ans[y];
		ans[y] = temp;
	}
	return ans;
}

bitset<BITSIZE> SA(bitset<BITSIZE> tempBit, MinCostFlowSolution * mcf) {

	if(nodesNum < 200) {
		may_server_num = 35;
	} else if (nodesNum < 500) {
		may_server_num = 60;
	} else {
		may_server_num = 120;
	}

	double t = T;
	bitset<BITSIZE> curPlan = tempBit;
	bitset<BITSIZE> newPlan;
	int P_L = 0;
	int P_F = 0;
	bitset<BITSIZE>  plan;
	long long cost = 1000000000;
	while(1) {

	    iterations++;
		long long curCost = 0;
		long long newCost = 0;
		for(int i = 0; i < ILOOP; i++) {
			newPlan = GetNext(curPlan);
			pair<pair<long long,long long>, bitset<BITSIZE> > pair1 = mcf->CalCost(curPlan, 0, false);
			curCost = pair1.first.first;
			pair1 = mcf->CalCost(newPlan, 0, false);
			newCost = pair1.first.first;
			double dE = newCost - curCost;

			if(dE < 0) {
				curPlan = newPlan;
				curCost = newCost;
				P_L = 0;
				P_F = 0;

				if(curCost < cost) {
				  plan = curPlan;
				  cost = curCost;
				  cout << "Cost: " << cost << ", Num: " << plan.count() << endl;
				}
			} else {
				double rd = get_random_real(0, 1);

				if(exp(-dE / t) > rd) {
					curPlan = newPlan;
					curCost = newCost;
				}
				P_L++;
			}
			if(P_L > LIMIT) {
				P_F++;
				break;
			}
		}

		if(P_F > OLOOP || iterations > Iterations){
			puts("bbbb");
			break;
		}
		//t *= DELTA;
		t = T * 1.0 /(1 + iterations);
	}
	cout<<"iterations :"<<iterations<<endl;
	return plan;
}

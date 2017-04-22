#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define TIME_LIMIT 85
#define EPS   1e-8
#define DELTA 0.98
// 20
#define T     50

#define LIMIT 200
#define OLOOP 18200
#define ILOOP 50 // 内层
#define Iterations 22000

#define P_GONE 0.8
#define P_SWAP 0.8
#define P_FLIP 0.7
#define P_CROSS 0.7

int may_server_num;

using namespace std;

int iterations = 0;
Individual GetNext(Individual ans) {

	if(int(ans.count()) > may_server_num) {
		double pGone = get_random_real(0, 1);
		if (pGone < P_GONE) {
			int postion = get_random_int(0, nodesNum - 1);
			if(ans.gen[postion]){ // 为 1 则反转
				ans.gen[postion] = 0;
			}
		}
	} else {
		double pFlip = get_random_real(0, 1);
		if(pFlip < P_FLIP) {
			int postion = get_random_int(0, nodesNum - 1);
			int le = get_random_int(1, serverLevel);
			if (ans.gen[postion])
				ans.gen[postion] = 0;
			else
				ans.gen[postion] = le;
		}
	}
	double pSwap = get_random_real(0, 1);
	if(pSwap < P_SWAP) {
		int x = get_random_int(0, nodesNum - 1);
		//int y = get_random_int(0, nodesNum - 1);
		int y;
		if (x == 0) {
			y = x + 1;
		} else {
			y = x - 1;
		}
		int temp = ans.gen[x];
		ans.gen[x] = ans.gen[y];
		ans.gen[y] = temp;
	}
	/*
	double pCross = get_random_real(0, 1);
	if (pCross < P_CROSS) {
		int x = get_random_int(0, nodesNum - 1);
		bitset<BITSIZE> bbb;
		bbb.reset();
		for (int i = 0; i < nodesNum - x; i++) {
			bbb[i] = ans[x+i];
		}
		for (int i = nodesNum - x; i < nodesNum; i++) {
			bbb[i] = ans[i - nodesNum + x];
		}
		ans = bbb;

	}
	*/
	return ans;
}

Individual SA(Individual tempIn, MinCostFlowSolution * mcf) {

	// int startSA = clock();

	if(nodesNum < 800) {
		may_server_num = 65;
	} else {
		may_server_num = 120;
	}

	double t = T;
	Individual curPlan = tempIn;
	Individual newPlan;
	int P_L = 0;
	int P_F = 0;
	Individual  plan;
	long long cost = 1000000000;
	while(iterations < Iterations) {

	    iterations++;
		long long curCost = 0;
		long long newCost = 0;
		for(int i = 0; i < ILOOP; i++) {
			newPlan = GetNext(curPlan);
			mcf->CalCost(curPlan, 0, false);
			curCost = curPlan.cost;
			mcf->CalCost(newPlan, 0, false);
			newCost = newPlan.cost;
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
			if (double(clock() - start) / CLOCKS_PER_SEC > TIME_LIMIT) {
				return plan;
			}
		}

		if(P_F > OLOOP){
			printf("OLOOP: P_F = %d\n", P_F);
			break;
		}
		// t *= DELTA;
		// t = T * 1.0 /(1 + iterations);
	}
	cout<<"iterations :"<<iterations<<endl;
	return plan;
}

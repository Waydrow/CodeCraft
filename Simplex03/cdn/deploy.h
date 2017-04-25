#ifndef __ROUTE_H__
#define __ROUTE_H__

#include "lib_io.h"
#include <bitset>

using namespace std;

#define BITSIZE 3700


void deploy_server(char * graph[MAX_EDGE_NUM], int edge_num, char * filename);
int calCost(bitset<BITSIZE> gene, int mustCal, bool);


#endif

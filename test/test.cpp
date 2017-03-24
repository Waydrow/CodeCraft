#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <bitset>

using namespace std;

int g_is_first = 1;
int uniform_int(int a, int b) {
    if (g_is_first) {
        g_is_first = 0;
        srand((unsigned int)time(NULL));
    }

    return (int)((double)rand() / ((RAND_MAX + 1.0) / (b - a + 1.0)) + a);
}

void test(int a, int b =1) {
    cout<<b;
}
int main() {
    double a, b;
    cin >> a >> b;
    cout <<min(a, b,1);
}

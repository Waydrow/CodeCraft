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
    bitset<10> a,b,c;
    a.reset();
    a.set(8);
    a.set(9);
    for (int i = 0; i < 10; i++) {
        cout << a[i];
    }
    cout <<endl;
    b = a<<1;
    for (int i = 0; i < 10; i++) {
        cout << b[i];
    }
    cout <<endl;
    c = a^b;
    for (int i = 0; i < 10; i++) {
        cout <<c[i];
    }
    cout <<endl;
}

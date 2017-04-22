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
void array_func(int a[]) {
    a[0] = 111;
}
int main() {
    int a[10];
    for (int i = 0; i < 10; i++) {
        a[i] = i;
    }

    array_func(a);
    for (int i = 0; i < 10; i++) {
        cout << a[i] << " ";
    }
    cout << endl;
    int n;
    cin >> n;
    int b[n];
    b[1] = 1;
    cout << b[1];
}

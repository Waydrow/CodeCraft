#include<iostream>
using namespace std;

void printbinary(const unsigned int val)
{
    for(int i = 2; i >= 0; i--)
    {
        if(val & (1 << i))
            cout << "1";
        else
            cout << "0";
    }
    cout << endl;
}

int main()
{
    while(true) {
        int x;
        cin >> x;
        printbinary(x);
    }
    return 0;
}

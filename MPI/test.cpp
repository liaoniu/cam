#include <iostream>
using namespace std;


void setZero(int *a)
{
    *a = 0;
}


int main()
{
    int x = 1;
    setZero(&x);
    cout << x << endl;

    return 0;
}
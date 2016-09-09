using namespace std;
#include "time.h"
int time()
{
    clock_t start, finish;
    start = clock();
    //code
    finish = clock();
    ((finish-start)/CLOCKS_PER_SEC);
}

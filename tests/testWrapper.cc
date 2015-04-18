/*
	program to test parameters.cc and parameters.h
*/

#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;

int main() {

for( int k=0; k<5; k++) {
	int N = 2;
	cout << " waiting for " << N << " units of time" << endl;
	sleep(N);
}

return 0;
}

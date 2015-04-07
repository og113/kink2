/*-------------------------------------------------------------------------------------------------------------------------
 	quick test for whether or not icc writes iterators for vectors of custom classes
 -------------------------------------------------------------------------------------------------------------------------*/
 
#include <iostream>
#include <vector>

using namespace std;

class testClass;

ostream& operator<<(ostream& os, const testClass& tc);

class testClass {
public:
	testClass(): number(7) {}
	~testClass() {}
	int value() const { return number;}
	
	friend ostream& operator<<(ostream& os, const testClass& tc);
private:
	int number;
};

ostream& operator<<(ostream& os, const testClass& tc) {
	os << tc.value();
	return os;
}

int main() {
cout << "test whether icc compiler writes iterators for std::vectors of custom classes" << endl;

vector<testClass> v(5);

for(unsigned int j=0; j<v.size(); j++) {
	cout << v[j] << endl;
}

return 0;
}

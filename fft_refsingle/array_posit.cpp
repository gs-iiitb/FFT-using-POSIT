#include <iostream>
#include <vector>
#include <iostream>
#include <typeinfo>
#include <random>
#include <limits>
#include "/home/gps/Downloads/My_test_POSIT/posit/posit.hpp"
#include "/home/gps/Downloads/My_test_POSIT/tests/test_helpers.hpp"
#include "/home/gps/Downloads/My_test_POSIT/tests/posit_math_helpers.hpp"

template<size_t nbits, size_t es>
void GenerateArray(){
sw::unum::posit<nbits, es> pa;
sw::unum::posit<nbits, es> pb[20];
int i;
for (i=0;i<20;i++){
pb[i]=pa.set_raw_bits(i);
std::cout<<sw::unum::to_hex(pb[i].get())<<std::endl;
}
}

int main(){
std::cout<<"Can I store POSIT into an array"<<std::endl;

GenerateArray<4,0>();

return 0;
}


#include <crab/config.h>
#include <crab/common/types.hpp>
#include <crab/common/wrapint.hpp>

#include <iostream>

using namespace std;

int main (int argc, char *argv[]) {

  using namespace crab;

  {
    wrapint n1(5, 3);
    wrapint n2(7, 3);
    wrapint n6(1450, 3);
    wrapint n3 = n1 + n2;
    wrapint n4 = n1 - n2;
    wrapint n5 = n1 * n2;  

    std::cout << "Bitwidth=3\n";
    crab::outs () << n1 << "+" << n2 << "=" << n3 << "\n";
    crab::outs () << n1 << "-" << n2 << "=" << n4 << "\n";
    crab::outs () << n1 << "*" << n2 << "=" << n5 << "\n";
    crab::outs () << "1450 converted to " << n6 << "\n";
  }

  {
    wrapint n1(5, 8);
    wrapint n2(7, 8);
    
    wrapint n3 = n1 + n2;
    wrapint n4 = n1 - n2;
    wrapint n5 = n1 * n2;
    wrapint n6(224, 8);
    wrapint n7(2, 8);
    wrapint n8 = n6 / n7;
    wrapint n9 = n1 % n2;

    std::cout << "Bitwidth=8\n";
    crab::outs () << n1 << "+" << n2 << "=" << n3 << "\n";
    crab::outs () << n1 << "-" << n2 << "=" << n4 << "\n";
    crab::outs () << n1 << "*" << n2 << "=" << n5 << "\n";
    crab::outs () << n6 << "/" << n7 << "=" << n8 << "\n";
    crab::outs () << n1 << "%" << n2 << "=" << n9 << "\n";        
  }

  {
    wrapint n1(5, 64);
    wrapint n2(7, 64);
    
    wrapint n3 = n1 + n2;
    wrapint n4 = n1 - n2;
    wrapint n5 = n1 * n2;  

    std::cout << "Bitwidth=64\n";
    crab::outs () << n1 << "+" << n2 << "=" << n3 << "\n";
    crab::outs () << n1 << "-" << n2 << "=" << n4 << "\n";
    crab::outs () << n1 << "*" << n2 << "=" << n5 << "\n";
  }

  {
    wrapint n1(7, 4);
    wrapint n2(12, 4);
    
    std::cout << "n1=" << std::bitset<4>(n1.get_uint64_t()) << "\n";
    std::cout << "n2=" << std::bitset<4>(n2.get_uint64_t()) << "\n";
    
    wrapint n3 = n1.sext(4);
    wrapint n4 = n2.sext(4);
    wrapint n5 = n1.zext(4);
    wrapint n6 = n2.zext(4);

    std::cout << "sext of n1 to 8 bits " << std::bitset<8>(n3.get_uint64_t()) << "\n";
    std::cout << "sext of n2 to 8 bits " << std::bitset<8>(n4.get_uint64_t()) << "\n";
    std::cout << "zext of n1 to 8 bits " << std::bitset<8>(n5.get_uint64_t()) << "\n";
    std::cout << "zext of n2 to 8 bits " << std::bitset<8>(n6.get_uint64_t()) << "\n";
  }

  {
    // 1010 1100
    wrapint n(-84, 8);
    std::cout << "n=" << std::bitset<8>(n.get_uint64_t()) << "\n";

    // 0000 0000 1010 1100 
    wrapint n1 = n.zext(8);
    // 1111 1111 1010 1100     
    wrapint n2 = n.sext(8);

    std::cout << "zext of n to 16 bits " << std::bitset<16>(n1.get_uint64_t()) << "\n";
    std::cout << "sext of n to 16 bits " << std::bitset<16>(n2.get_uint64_t()) << "\n";
  }
  
  {
    wrapint n1(7, 32);
    wrapint n2(-999999, 32);
    
    std::cout << "n1=" << std::bitset<32>(n1.get_uint64_t()) << "\n";
    std::cout << "n2=" << std::bitset<32>(n2.get_uint64_t()) << "\n";
    
    wrapint n3 = n1.sext(31);
    wrapint n4 = n2.sext(31);
    wrapint n5 = n1.zext(31);
    wrapint n6 = n2.zext(31);
    
    std::cout << "sext of n1 to 63 bits " << std::bitset<63>(n3.get_uint64_t()) << "\n";
    std::cout << "sext of n2 to 63 bits " << std::bitset<63>(n4.get_uint64_t()) << "\n";
    std::cout << "zext of n1 to 63 bits " << std::bitset<63>(n5.get_uint64_t()) << "\n";
    std::cout << "zext of n2 to 63 bits " << std::bitset<63>(n6.get_uint64_t()) << "\n";
  }

  {
    wrapint n1(7, 32);
    wrapint n2(-999999, 32);
    
    std::cout << "n1=" << std::bitset<32>(n1.get_uint64_t()) << "\n";
    std::cout << "n2=" << std::bitset<32>(n2.get_uint64_t()) << "\n";
    
    wrapint n3 = n1.sext(32);
    wrapint n4 = n2.sext(32);
    wrapint n5 = n1.zext(32);
    wrapint n6 = n2.zext(32);

    std::cout << "sext of n1 to 64 bits " << std::bitset<64>(n3.get_uint64_t()) << "\n";
    std::cout << "sext of n2 to 64 bits " << std::bitset<64>(n4.get_uint64_t()) << "\n";
    std::cout << "zext of n1 to 64 bits " << std::bitset<64>(n5.get_uint64_t()) << "\n";
    std::cout << "zext of n2 to 64 bits " << std::bitset<64>(n6.get_uint64_t()) << "\n";
  }
  
  
  return 0;
}

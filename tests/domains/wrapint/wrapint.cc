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
  
  
  return 0;
}

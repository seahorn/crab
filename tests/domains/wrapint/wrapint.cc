#include <crab/config.h>
#include <crab/numbers/wrapint.hpp>

#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

  using namespace crab;

  {
    wrapint n1(5, 3);
    wrapint n2(7, 3);
    wrapint n6(1450, 3);
    wrapint n3 = n1 + n2;
    wrapint n4 = n1 - n2;
    wrapint n5 = n1 * n2;

    std::cout << "Bitwidth=3\n";
    crab::outs() << n1 << "+" << n2 << "=" << n3 << "\n";
    crab::outs() << n1 << "-" << n2 << "=" << n4 << "\n";
    crab::outs() << n1 << "*" << n2 << "=" << n5 << "\n";
    crab::outs() << "1450 converted to " << n6 << "\n";
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
    crab::outs() << n1 << "+" << n2 << "=" << n3 << "\n";
    crab::outs() << n1 << "-" << n2 << "=" << n4 << "\n";
    crab::outs() << n1 << "*" << n2 << "=" << n5 << "\n";
    crab::outs() << n6 << "/" << n7 << "=" << n8 << "\n";
    crab::outs() << n1 << "%" << n2 << "=" << n9 << "\n";
  }

  {
    wrapint n1(5, 64);
    wrapint n2(7, 64);

    wrapint n3 = n1 + n2;
    wrapint n4 = n1 - n2;
    wrapint n5 = n1 * n2;

    std::cout << "Bitwidth=64\n";
    crab::outs() << n1 << "+" << n2 << "=" << n3 << "\n";
    crab::outs() << n1 << "-" << n2 << "=" << n4 << "\n";
    crab::outs() << n1 << "*" << n2 << "=" << n5 << "\n";
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

    std::cout << "sext of n1 to 8 bits " << std::bitset<8>(n3.get_uint64_t())
              << "\n";
    std::cout << "sext of n2 to 8 bits " << std::bitset<8>(n4.get_uint64_t())
              << "\n";
    std::cout << "zext of n1 to 8 bits " << std::bitset<8>(n5.get_uint64_t())
              << "\n";
    std::cout << "zext of n2 to 8 bits " << std::bitset<8>(n6.get_uint64_t())
              << "\n";
  }

  {
    // 1010 1100
    wrapint n(-84, 8);
    std::cout << "n=" << std::bitset<8>(n.get_uint64_t()) << "\n";

    // 0000 0000 1010 1100
    wrapint n1 = n.zext(8);
    // 1111 1111 1010 1100
    wrapint n2 = n.sext(8);

    std::cout << "zext of n to 16 bits " << std::bitset<16>(n1.get_uint64_t())
              << "\n";
    std::cout << "sext of n to 16 bits " << std::bitset<16>(n2.get_uint64_t())
              << "\n";
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

    std::cout << "sext of n1 to 63 bits " << std::bitset<63>(n3.get_uint64_t())
              << "\n";
    std::cout << "sext of n2 to 63 bits " << std::bitset<63>(n4.get_uint64_t())
              << "\n";
    std::cout << "zext of n1 to 63 bits " << std::bitset<63>(n5.get_uint64_t())
              << "\n";
    std::cout << "zext of n2 to 63 bits " << std::bitset<63>(n6.get_uint64_t())
              << "\n";
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

    std::cout << "sext of n1 to 64 bits " << std::bitset<64>(n3.get_uint64_t())
              << "\n";
    std::cout << "sext of n2 to 64 bits " << std::bitset<64>(n4.get_uint64_t())
              << "\n";
    std::cout << "zext of n1 to 64 bits " << std::bitset<64>(n5.get_uint64_t())
              << "\n";
    std::cout << "zext of n2 to 64 bits " << std::bitset<64>(n6.get_uint64_t())
              << "\n";
  }

  {                    // shift operators
    wrapint n1(13, 4); //
    wrapint n2(2, 4);  //

    wrapint n3 = n1 << n2;
    std::cout << std::bitset<4>(n1.get_uint64_t()) << " << "
              << std::bitset<4>(n2.get_uint64_t()) << "="
              << std::bitset<4>(n3.get_uint64_t()) << "\n";
    wrapint n4 = n1.ashr(n2);
    std::cout << std::bitset<4>(n1.get_uint64_t()) << " >>_a "
              << std::bitset<4>(n2.get_uint64_t()) << "="
              << std::bitset<4>(n4.get_uint64_t()) << "\n";
    wrapint n5 = n1.lshr(n2);
    std::cout << std::bitset<4>(n1.get_uint64_t()) << " >>_l "
              << std::bitset<4>(n2.get_uint64_t()) << "="
              << std::bitset<4>(n5.get_uint64_t()) << "\n";
  }

  { // signed vs unsigned division
    wrapint n1(9, 4);
    wrapint n2(2, 4);
    wrapint n3 = n1.sdiv(n2);
    wrapint n4 = n1.udiv(n2);
    std::cout << std::bitset<4>(n1.get_uint64_t()) << " /_s "
              << std::bitset<4>(n2.get_uint64_t()) << " ="
              << std::bitset<4>(n3.get_uint64_t()) << "\n";
    std::cout << std::bitset<4>(n1.get_uint64_t()) << " /_u "
              << std::bitset<4>(n2.get_uint64_t()) << " ="
              << std::bitset<4>(n4.get_uint64_t()) << "\n";
  }

  {                    // signed remainder
    wrapint n1(7, 4);  // 7
    wrapint n2(3, 4);  // 3
    wrapint n3(13, 4); // -3
    wrapint n4(9, 4);  // -7

    wrapint n5 = n1.srem(n2);
    wrapint n6 = n1.srem(n3);
    wrapint n7 = n4.srem(n2);
    wrapint n8 = n4.srem(n3);

    // 7 %_s 3 = 1
    std::cout << std::bitset<4>(n1.get_uint64_t()) << " %_s "
              << std::bitset<4>(n2.get_uint64_t()) << " ="
              << std::bitset<4>(n5.get_uint64_t());
    crab::outs() << " (" << n5.get_signed_bignum() << ") \n";

    // 7 %_s -3 = 1
    std::cout << std::bitset<4>(n1.get_uint64_t()) << " %_s "
              << std::bitset<4>(n3.get_uint64_t()) << " ="
              << std::bitset<4>(n6.get_uint64_t());
    crab::outs() << " (" << n6.get_signed_bignum() << ") \n";

    // -7 %_s 3 = -1
    std::cout << std::bitset<4>(n4.get_uint64_t()) << " %_s "
              << std::bitset<4>(n2.get_uint64_t()) << " ="
              << std::bitset<4>(n7.get_uint64_t());
    crab::outs() << " (" << n7.get_signed_bignum() << ") \n";

    // -7 %_s -3 = -1
    std::cout << std::bitset<4>(n4.get_uint64_t()) << " %_s "
              << std::bitset<4>(n3.get_uint64_t()) << " ="
              << std::bitset<4>(n8.get_uint64_t());
    crab::outs() << " (" << n8.get_signed_bignum() << ") \n";
  }

  {                    // unsigned remainder
    wrapint n1(7, 4);  // 7
    wrapint n2(3, 4);  // 3
    wrapint n3(13, 4); // 13
    wrapint n4(9, 4);  // 9

    wrapint n5 = n1.urem(n2);
    wrapint n6 = n1.urem(n3);
    wrapint n7 = n4.urem(n2);
    wrapint n8 = n4.urem(n3);

    // 7 %_u 3 = 1
    std::cout << std::bitset<4>(n1.get_uint64_t()) << " %_u "
              << std::bitset<4>(n2.get_uint64_t()) << " ="
              << std::bitset<4>(n5.get_uint64_t()) << "\n";

    // 7 %_s 13 = 7
    std::cout << std::bitset<4>(n1.get_uint64_t()) << " %_u "
              << std::bitset<4>(n3.get_uint64_t()) << " ="
              << std::bitset<4>(n6.get_uint64_t()) << "\n";

    // 9 %_s 3 = 0
    std::cout << std::bitset<4>(n4.get_uint64_t()) << " %_s "
              << std::bitset<4>(n2.get_uint64_t()) << " ="
              << std::bitset<4>(n7.get_uint64_t()) << "\n";

    // 9 %_s 13 = 9
    std::cout << std::bitset<4>(n4.get_uint64_t()) << " %_s "
              << std::bitset<4>(n3.get_uint64_t()) << " ="
              << std::bitset<4>(n8.get_uint64_t()) << "\n";
  }

  {
    wrapint v1 = wrapint::get_unsigned_max(64);
    wrapint v2 = wrapint::get_unsigned_min(64);
    wrapint v3 = wrapint::get_signed_max(64);
    wrapint v4 = wrapint::get_signed_min(64);
    wrapint v5 = wrapint::get_unsigned_max(32);
    wrapint v6 = wrapint::get_unsigned_min(32);
    wrapint v7 = wrapint::get_signed_max(32);
    wrapint v8 = wrapint::get_signed_min(32);
    wrapint v9 = wrapint::get_unsigned_max(40);
    wrapint v10 = wrapint::get_unsigned_min(40);
    wrapint v11 = wrapint::get_signed_max(40);
    wrapint v12 = wrapint::get_signed_min(40);
    wrapint v13 = wrapint::get_unsigned_max(63);
    wrapint v14 = wrapint::get_unsigned_min(63);
    wrapint v15 = wrapint::get_signed_max(63);
    wrapint v16 = wrapint::get_signed_min(63);

    crab::outs() << "UMAX(64)=" << v1 << "\n";
    std::cout << std::bitset<64>(v1.get_uint64_t()) << "\n";
    crab::outs() << "UMIN(64)=" << v2 << "\n";
    std::cout << std::bitset<64>(v2.get_uint64_t()) << "\n";
    crab::outs() << "SMAX(64)=" << v3 << "\n";
    std::cout << std::bitset<64>(v3.get_uint64_t()) << "\n";
    crab::outs() << "SMIN(64)=" << v4 << "\n";
    std::cout << std::bitset<64>(v4.get_uint64_t()) << "\n";
    crab::outs() << "UMAX(32)=" << v5 << "\n";
    std::cout << std::bitset<32>(v5.get_uint64_t()) << "\n";
    crab::outs() << "UMIN(32)=" << v6 << "\n";
    std::cout << std::bitset<32>(v6.get_uint64_t()) << "\n";
    crab::outs() << "SMAX(32)=" << v7 << "\n";
    std::cout << std::bitset<32>(v7.get_uint64_t()) << "\n";
    crab::outs() << "SMIN(32)=" << v8 << "\n";
    std::cout << std::bitset<32>(v8.get_uint64_t()) << "\n";
    crab::outs() << "UMAX(40)=" << v9 << "\n";
    std::cout << std::bitset<40>(v9.get_uint64_t()) << "\n";
    crab::outs() << "UMIN(40)=" << v10 << "\n";
    std::cout << std::bitset<40>(v10.get_uint64_t()) << "\n";
    crab::outs() << "SMAX(40)=" << v11 << "\n";
    std::cout << std::bitset<40>(v11.get_uint64_t()) << "\n";
    crab::outs() << "SMIN(40)=" << v12 << "\n";
    std::cout << std::bitset<40>(v12.get_uint64_t()) << "\n";
    crab::outs() << "UMAX(63)=" << v13 << "\n";
    std::cout << std::bitset<63>(v13.get_uint64_t()) << "\n";
    crab::outs() << "UMIN(63)=" << v14 << "\n";
    std::cout << std::bitset<63>(v14.get_uint64_t()) << "\n";
    crab::outs() << "SMAX(63)=" << v15 << "\n";
    std::cout << std::bitset<63>(v15.get_uint64_t()) << "\n";
    crab::outs() << "SMIN(63)=" << v16 << "\n";
    std::cout << std::bitset<63>(v16.get_uint64_t()) << "\n";
  }

  return 0;
}

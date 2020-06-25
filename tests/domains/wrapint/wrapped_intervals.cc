#include "../../program_options.hpp"
#include <crab/config.h>
#include <crab/domains/wrapped_interval_domain.hpp>
#include <crab/numbers/wrapint.hpp>

using namespace std;
using namespace crab::domains;

int main(int argc, char *argv[]) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  using wrapped_interval_t = wrapped_interval<ikos::z_number>;

  {
    wrapped_interval_t i1(crab::wrapint(-10, 8), crab::wrapint(10, 8));
    wrapped_interval_t i2(crab::wrapint(0, 8), crab::wrapint(20, 8));

    crab::outs() << "i1=" << i1 << "\n";
    crab::outs() << "i2=" << i2 << "\n";
    crab::outs() << "-i2=" << (-i2) << "\n";
    crab::outs() << "i1 & i2: " << (i1 & i2) << "\n";
    crab::outs() << "i1 | i2: " << (i1 | i2) << "\n";

    wrapped_interval_t i3(crab::wrapint(0, 8), crab::wrapint(-10, 8));
    wrapped_interval_t i4(crab::wrapint(-15, 8), crab::wrapint(5, 8));
    crab::outs() << i3 << " & " << i4 << ": " << (i3 & i4) << "\n";
    crab::outs() << i3 << " | " << i4 << ": " << (i3 | i4) << "\n";

    wrapped_interval_t urk(crab::wrapint(1642571628, 32),
                           crab::wrapint(700177772, 32));
    crab::outs() << "urk: " << urk << "\n";
    crab::outs() << "urk.is_top(): " << urk.is_top() << "\n";
    crab::outs() << "urk == Top: " << (urk == wrapped_interval_t::top())
                 << "\n";
  }

  {
    // mutiplication
    wrapped_interval_t i1(crab::wrapint(15, 4), crab::wrapint(9, 4));
    wrapped_interval_t i2(crab::wrapint(0, 4), crab::wrapint(1, 4));
    wrapped_interval_t i3 = i1 * i2;
    crab::outs() << i1 << " * " << i2 << " = " << i3 << "\n";
  }

  {
    // signed and unsigned division
    // [4,7] /_s  [-2,3] = [1,-2]
    // [4,7] /_u  [-2,3] = [0,7]
    //    [4,7] /_u [-2,-1] U [4,7] /_u [1,3] =
    //    [4,7] /_u [14,15] U [4,7] /_u [1,3] =
    //    [4/15, 7/14] U [4/3, 7/1] =
    //    [0,0] U [1,7] = [0,7]
    wrapped_interval_t i1(crab::wrapint(4, 4), crab::wrapint(7, 4));
    wrapped_interval_t i2(crab::wrapint(14, 4), crab::wrapint(3, 4));
    crab::outs() << i1 << " /_s " << i2 << " = " << i1.SDiv(i2) << "\n";
    crab::outs() << i1 << " /_u " << i2 << " = " << i1.UDiv(i2) << "\n";
    wrapped_interval_t i5(crab::wrapint(4, 4), crab::wrapint(7, 4));
    wrapped_interval_t i6(crab::wrapint(2, 4), crab::wrapint(2, 4));
    crab::outs() << i5 << " /_s " << i6 << " = " << i5.SDiv(i6) << "\n";
    crab::outs() << i5 << " /_u " << i6 << " = " << i5.UDiv(i6) << "\n";
  }

  return 0;
}

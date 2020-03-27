#include "../program_options.hpp"
#include "../common.hpp"

using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

int main (int argc, char** argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }

  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var v(vfac["v"], crab::INT_TYPE, 32);    

  {
    crab::outs() << " ===== " << z_interval_domain_t::getDomainName() << " ====\n";
    z_interval_domain_t inv;
    inv += z_lin_cst_t(x <= 0);
    inv += z_lin_cst_t(y == 5);
    crab::outs() << "Before renaming " << inv << "\n";

    {
      z_interval_domain_t tmp(inv);
      tmp.rename({x,y}, {w,z});
      crab::outs() << "After rename x with w and y with z=" << tmp << "\n";
    }
    {
      z_interval_domain_t tmp(inv);
      tmp.rename({x}, {w});
      crab::outs() << "After rename x with w=" << tmp << "\n";
    }
    // {
    //   z_interval_domain_t tmp(inv);
    //   tmp.rename({x}, {y});
    //   crab::outs() << "After rename x with y=" << tmp << "\n";
    //   // This should raise an error:
    //   // CRAB ERROR: rename expects y to be unconstrained
    // }
  }

#ifdef HAVE_APRON  
  {
    crab::outs() << " ===== " << z_pk_apron_domain_t::getDomainName() << " ====\n";    
    z_pk_apron_domain_t inv;
    inv += z_lin_cst_t(x <= 0);
    inv += z_lin_cst_t(y == 5);
    crab::outs() << "Before renaming " << inv << "\n";

    {
      z_pk_apron_domain_t tmp(inv);
      tmp.rename({x,y}, {w,z});
      crab::outs() << "After rename x with w and y with z=" << tmp << "\n";
    }
    {
      z_pk_apron_domain_t tmp(inv);
      tmp.rename({x}, {w});
      crab::outs() << "After rename x with w=" << tmp << "\n";
    }
    // {
    //   z_pk_apron_domain_t tmp(inv);
    //   tmp.rename({x}, {y});
    //   crab::outs() << "After rename x with y=" << tmp << "\n";
    //   // This should raise an error:
    //   //  CRAB ERROR: rename expects y to be unconstrained
    // }
  }
#endif   
#ifdef HAVE_ELINA
  {
    crab::outs() << " ===== " << z_pk_elina_domain_t::getDomainName() << " ====\n";    
    z_pk_elina_domain_t inv;
    inv += z_lin_cst_t(x <= 0);
    inv += z_lin_cst_t(y == 5);
    crab::outs() << "Before renaming " << inv << "\n";

    {
      z_pk_elina_domain_t tmp(inv);
      tmp.rename({x,y}, {w,z});
      crab::outs() << "After rename x with w and y with z=" << tmp << "\n";
    }
    {
      z_pk_elina_domain_t tmp(inv);
      tmp.rename({x}, {w});
      crab::outs() << "After rename x with w=" << tmp << "\n";
    }
    // {
    //   z_pk_elina_domain_t tmp(inv);
    //   tmp.rename({x}, {y});
    //   crab::outs() << "After rename x with y=" << tmp << "\n";
    //   // This should raise an error:
    //   //  CRAB ERROR: rename expects y to be unconstrained
    // }
  }
#endif   

  {
    crab::outs() << " ===== " << z_sdbm_domain_t::getDomainName() << " ====\n";        
    z_sdbm_domain_t inv;
    inv += z_lin_cst_t(x <= 0);
    inv += z_lin_cst_t(y == 5);
    crab::outs() << "Before renaming " << inv << "\n";

    {
      z_sdbm_domain_t tmp(inv);
      tmp.rename({x,y}, {w,z});
      crab::outs() << "After rename x with w and y with z=" << tmp << "\n";
    }
    {
      z_sdbm_domain_t tmp(inv);
      tmp.rename({x}, {w});
      crab::outs() << "After rename x with w=" << tmp << "\n";
    }
    // {
    //   z_sdbm_domain_t tmp(inv);
    //   tmp.rename({x}, {y});
    //   crab::outs() << "After rename x with y=" << tmp << "\n";
    //   // This should raise an error:
    //   //  CRAB ERROR: rename expects y to be unconstrained
    // }
  }
  
  {
    crab::outs() << " ===== " << z_dbm_domain_t::getDomainName() << " ====\n";            
    z_dbm_domain_t inv;
    inv += z_lin_cst_t(x <= 0);
    inv += z_lin_cst_t(y == 5);
    crab::outs() << "Before renaming " << inv << "\n";

    {
      z_dbm_domain_t tmp(inv);
      tmp.rename({x,y}, {w,z});
      crab::outs() << "After rename x with w and y with z=" << tmp << "\n";
    }
    {
      z_dbm_domain_t tmp(inv);
      tmp.rename({x}, {w});
      crab::outs() << "After rename x with w=" << tmp << "\n";
    }
    // {
    //   z_dbm_domain_t tmp(inv);
    //   tmp.rename({x}, {y});
    //   crab::outs() << "After rename x with y=" << tmp << "\n";
    //   // This should raise an error:
    //   //  CRAB ERROR: rename expects y to be unconstrained    
    // }
  }

  {
    crab::outs() << " ===== " << z_term_domain_t::getDomainName() << " ====\n";            
    z_term_domain_t inv;
    inv += z_lin_cst_t(x <= 0);
    inv += z_lin_cst_t(y == 5);
    crab::outs() << "Before renaming " << inv << "\n";

    {
      z_term_domain_t tmp(inv);
      tmp.rename({x,y}, {w,z});
      crab::outs() << "After rename x with w and y with z=" << tmp << "\n";
    }
    {
      z_term_domain_t tmp(inv);
      tmp.rename({x}, {w});
      crab::outs() << "After rename x with w=" << tmp << "\n";
    }
    // {
    //   z_term_domain_t tmp(inv);
    //   tmp.rename({x}, {y});
    //   crab::outs() << "After rename x with y=" << tmp << "\n";
    //   // This should raise an error:
    //   //  CRAB ERROR: rename expects y to be unconstrained    
    // }
  }
  
  return 0;
}

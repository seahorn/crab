#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
cfg_t* prog1 (variable_factory_t &vfac)  {

  // Definining program variables
  z_var i (vfac ["i"]);
  z_var k (vfac ["k"]);
  z_var nd (vfac ["nd"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry","ret");
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& bb1   = cfg->insert ("bb1");
  basic_block_t& bb1_t = cfg->insert ("bb1_t");
  basic_block_t& bb1_f = cfg->insert ("bb1_f");
  basic_block_t& bb2   = cfg->insert ("bb2");
  basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign (k, 0);
  entry.assign (i, 0);
  bb1_t.assume (i <= 99);
  bb1_f.assume (i >= 100);
  bb2.add(i, i, 1);
  bb2.add(k, k, 1);
  return cfg;
}

cfg_t* prog2 (variable_factory_t &vfac) 
{

  cfg_t* cfg = new cfg_t("loop1_entry","ret"); 
  basic_block_t& loop1_entry = cfg->insert ("loop1_entry");
  basic_block_t& loop1_bb1   = cfg->insert ("loop1_bb1");
  basic_block_t& loop1_bb1_t = cfg->insert ("loop1_bb1_t");
  basic_block_t& loop1_bb1_f = cfg->insert ("loop1_bb1_f");
  basic_block_t& loop1_bb2   = cfg->insert ("loop1_bb2");
  basic_block_t& loop2_entry = cfg->insert ("loop2_entry");
  basic_block_t& loop2_bb1   = cfg->insert ("loop2_bb1");
  basic_block_t& loop2_bb1_t = cfg->insert ("loop2_bb1_t");
  basic_block_t& loop2_bb1_f = cfg->insert ("loop2_bb1_f");
  basic_block_t& loop2_bb2   = cfg->insert ("loop2_bb2");
  basic_block_t& ret         = cfg->insert ("ret");

  loop1_entry >> loop1_bb1;
  loop1_bb1 >> loop1_bb1_t; loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >> loop1_bb2; loop1_bb2 >> loop1_bb1; loop1_bb1_f >> loop2_entry;

  loop2_entry >> loop2_bb1;
  loop2_bb1 >> loop2_bb1_t; loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb2; loop2_bb2 >> loop2_bb1; loop2_bb1_f >> ret;

  z_var i(vfac["i"]);
  z_var j(vfac["j"]);
  z_var k(vfac["k"]);

  loop1_entry.assign (i, 0);
  loop1_entry.assign (k, 30);
  loop1_bb1_t.assume (i <= 9);
  loop1_bb1_f.assume (i >= 10);
  loop1_bb2.add (i, i, 1);

  loop2_entry.assign (j, 0);
  loop2_bb1_t.assume (j <= 9);
  loop2_bb1_f.assume (j >= 10);
  loop2_bb2.add (j, j, 1);
  return cfg;
}

cfg_t* prog3 (variable_factory_t &vfac) 
{

  cfg_t* cfg = new cfg_t("entry","ret");
  basic_block_t& entry       = cfg->insert ("entry");
  basic_block_t& loop1_head  = cfg->insert ("loop1_head");
  basic_block_t& loop1_t     = cfg->insert ("loop1_t");
  basic_block_t& loop1_f     = cfg->insert ("loop1_f");
  basic_block_t& loop1_body  = cfg->insert ("loop1_body");

  basic_block_t& loop1_body_t  = cfg->insert ("loop1_body_t");
  basic_block_t& loop1_body_f  = cfg->insert ("loop1_body_f");
  basic_block_t& loop1_body_x  = cfg->insert ("loop1_body_x");

  basic_block_t& cont        = cfg->insert ("cont");
  basic_block_t& loop2_head  = cfg->insert ("loop2_head");
  basic_block_t& loop2_t     = cfg->insert ("loop2_t");
  basic_block_t& loop2_f     = cfg->insert ("loop2_f");
  basic_block_t& loop2_body  = cfg->insert ("loop2_body");
  basic_block_t& ret         = cfg->insert ("ret");

  entry >> loop1_head;
  loop1_head >> loop1_t; 
  loop1_head >> loop1_f; 
  loop1_t >>    loop1_body; 

  loop1_body >> loop1_body_t;
  loop1_body >> loop1_body_f;
  loop1_body_t >> loop1_body_x;
  loop1_body_f >> loop1_body_x;
  loop1_body_x >> loop1_head;

  loop1_f >> cont;
  cont >> loop2_head;
  loop2_head >> loop2_t; 
  loop2_head >> loop2_f; 
  loop2_t >>    loop2_body; 
  loop2_body >> loop2_head;
  loop2_f >> ret;
  
  z_var i(vfac["i"]);

  entry.assign (i, 0);
  loop1_t.assume (i <= 10);
  loop1_f.assume (i >= 11);
  loop1_body.add (i, i, 1);

  loop1_body_t.assume (i >= 9);
  loop1_body_t.assign (i , 0);
  loop1_body_f.assume (i <= 8);

  loop2_t.assume (i <= 100);
  loop2_f.assume (i >= 101);
  loop2_body.sub (i, i, 1);
  return cfg;
}

cfg_t* prog4 (variable_factory_t &vfac) 
{

  cfg_t* cfg = new cfg_t("entry","ret");
  basic_block_t& entry      = cfg->insert ("entry");
  basic_block_t& loop_head  = cfg->insert ("loop_head");
  basic_block_t& loop_t     = cfg->insert ("loop_t");
  basic_block_t& loop_f     = cfg->insert ("loop_f");
  basic_block_t& loop_body  = cfg->insert ("loop_body");
  basic_block_t& ret        = cfg->insert ("ret");

  entry >> loop_head;
  loop_head >> loop_t; 
  loop_head >> loop_f; 
  loop_t >> loop_body; 
  loop_body >> loop_head;
  loop_f >> ret;

  z_var i(vfac["i"]);
  z_var p(vfac["p"]);

  entry.assign (i, 0);
  entry.assign (p, 0);

  loop_t.assume (i <= 9);
  loop_f.assume (i >= 10);
  loop_body.add (i, i, 1);
  loop_body.add (p, p, 4);

  return cfg;
}

/* Example of how to build a CFG */
cfg_t* prog5 (variable_factory_t &vfac)  {

  // Definining program variables
  z_var i (vfac ["i"]);
  z_var k (vfac ["k"]);
  z_var nd (vfac ["nd"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry","ret");
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& bb1   = cfg->insert ("bb1");
  basic_block_t& bb1_t = cfg->insert ("bb1_t");
  basic_block_t& bb1_f = cfg->insert ("bb1_f");
  basic_block_t& bb2   = cfg->insert ("bb2");
  basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign (k, 0);
  entry.assign (i, 0);
  bb1_t.assume (i != 9);
  bb1_f.assume (i == 9);
  bb2.add(i, i, 1);
  bb2.add(k, k, 1);
  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;

  cfg_t* cfg1 = prog1 (vfac);      
  cfg_t* cfg2 = prog2 (vfac);
  cfg_t* cfg3 = prog3 (vfac);
  cfg_t* cfg4 = prog4 (vfac);
  cfg_t* cfg5 = prog5 (vfac);

  cfg1->simplify (); // this is optional
  crab::outs() << *cfg1 << "\n";
  run<boxes_domain_t>(cfg1, vfac, false, 10,2,20);

  cfg2->simplify (); // this is optional
  crab::outs() << *cfg2 << "\n";
  run<boxes_domain_t>(cfg2, vfac, false, 10,2,20);
  
  cfg3->simplify (); // this is optional
  crab::outs() << *cfg3 << "\n";
  run<boxes_domain_t>(cfg3, vfac, false, 10,2,20);

  cfg4->simplify (); // this is optional
  crab::outs() << *cfg4 << "\n";
  run<boxes_domain_t>(cfg4, vfac, false, 10,2,20);

  crab::outs() << *cfg5 << "\n";
  run<boxes_domain_t>(cfg5, vfac, false, 10,2,20);

  delete cfg1;
  delete cfg2;
  delete cfg3;
  delete cfg4;
  delete cfg5;

  { 
    crab::outs() << "Testing some boxes operations ...\n";
    varname_t x = vfac["x"];
    varname_t y = vfac["y"];
    varname_t z = vfac["z"];

    boxes_domain_t inv1 = boxes_domain_t::top ();

    inv1.assign (y, 6);
    inv1.assign (z, 7);

    boxes_domain_t inv2 = boxes_domain_t::top ();
    inv2.assign (y, 3);
    inv2.assign (z, 4);

    boxes_domain_t inv3 = inv1 | inv2;

    crab::outs() << inv3 << "\n";
   
    crab::outs() << x << ":=" << y << " + " << z << "= \n";
    inv3.apply (OP_ADDITION, x, y, z);
    crab::outs() << inv3 << "\n";

    crab::outs() << x << ":=" << y << " - " << z << "= \n";
    inv3.apply (OP_SUBTRACTION, x, y, z);
    crab::outs() << inv3 << "\n";

    crab::outs() << x << ":=" << y << " * " << z << "= \n";
    inv3.apply (OP_MULTIPLICATION, x, y, z);
    crab::outs() << inv3 << "\n";

    crab::outs() << x << ":=" << y << " / " << z << "= \n";
    inv3.apply (OP_DIVISION, x, y, z);
    crab::outs() << inv3 << "\n";

    boxes_domain_t inv4 = boxes_domain_t::top ();    
    z_var cx (vfac["x"]);
    z_var cy (vfac["y"]);
    z_var cz (vfac["z"]);

    inv4 +=  (cx >= cy);
    crab::outs() << "Added x >= y \n" << inv4 << "\n";    

    inv4 +=  (cx != 9);
    crab::outs() << "Added x != 9\n" << inv4 << "\n";    

    inv4 +=  (cy >=  9);
    crab::outs() << "Added y >= 9\n" << inv4 << "\n";    

    inv4 +=  (cy <=  9);
    crab::outs() << "Added y <= 9\n" << inv4 << "\n";    

    inv4 +=  (cz >= 10); 
    crab::outs() << "Added z > 9\n" << inv4 << "\n";    

    inv4 +=  (cz <= 8);
    crab::outs() << "Added z < 9\n" << inv4 << "\n";    
  }
  return 0;
}

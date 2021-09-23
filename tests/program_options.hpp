#pragma once

#include <boost/program_options.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/domains/abstract_domain_params.hpp>

#include <iostream>

namespace crab_tests {

int parse_user_options(int argc, char **argv, bool &stats_enabled) {
  boost::program_options::options_description po("Test Options");
  po.add_options()("help", "Print help message and exit");
  po.add_options()("log",
                   boost::program_options::value<std::vector<std::string>>(),
                   "Enable specified log level");
  po.add_options()("domain-param",
                   boost::program_options::value<std::vector<std::string>>(),
                   "Set abstract domain parameter: arg must be \"param=val\"");
  po.add_options()("display-domain-params",
		   "Display abstract domain parameter values");
  po.add_options()("verbose", boost::program_options::value<unsigned>(),
                   "Enable verbosity level");
  po.add_options()("stats", boost::program_options::bool_switch(&stats_enabled),
                   "Enable stats");
  po.add_options()("disable-warnings", "Disable warning messages");
  po.add_options()("sanity", "Enable sanity checks");
  boost::program_options::options_description cmmdline_options;
  cmmdline_options.add(po);
  boost::program_options::variables_map vm;
  boost::program_options::positional_options_description p;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv)
          .options(cmmdline_options)
          .positional(p)
          .run(),
      vm);
  boost::program_options::notify(vm);
  if (vm.count("help")) {
    std::cout << po << "n";
    return 0;
  }
  if (vm.count("log")) {
    std::vector<std::string> loggers = vm["log"].as<std::vector<std::string>>();
    for (unsigned int i = 0; i < loggers.size(); i++)
      crab::CrabEnableLog(loggers[i]);
  }
  
  if (vm.count("domain-param")) {
    std::vector<std::string> parameters = vm["domain-param"].as<std::vector<std::string>>();
    std::string delimiter("=");
    for (unsigned int i = 0; i < parameters.size(); i++) {
      std::string s = parameters[i];
      int pos = s.find(delimiter);
      std::string param = s.substr(0, pos);
      std::string val = s.substr(pos+delimiter.length());
      crab::domains::crab_domain_params_man::get().set_param(param, val);
    }
  }
  // print after set all domain-param
  if (vm.count("display-domain-params")) {
    crab::domains::crab_domain_params_man::get().write(crab::outs());
  }
  
  if (vm.count("verbose")) {
    crab::CrabEnableVerbosity(vm["verbose"].as<unsigned>());
  }
  if (vm.count("disable-warnings")) {
    crab::CrabEnableWarningMsg(false);
  }
  if (vm.count("sanity")) {
    crab::CrabEnableSanityChecks(true);
  }
  crab::CrabEnableStats();
  return 1;
}

} // namespace crab_tests

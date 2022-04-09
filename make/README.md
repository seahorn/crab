# How to use Crab via Makefile #

This directory contains a sample `Makefile` and two programs:

- `domain.cpp`: use Crab C++ API to create an **abstract domain
  object** and perform some abstract operations.

- `analysis.cpp`: use Crab C++ API to create a Control Flow Graph
  (CFG) and an **analysis object**, and perform analysis on the CFG.

There are also two subdirectories `apron` and `elina` that show how to
use Apron and Elina libraries. Currently, they are not compatible with
each other so that's why they must be built separately.

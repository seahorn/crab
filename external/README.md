# How to use Crab from an external application #

This directory contains a simple `Makefile` and two minimalistic
programs:

- `domain.cpp`: use Crab C++ API to create an **abstract domain
  object** and perform some abstract operations.

- `analysis.cpp`: use Crab C++ API to create a Control Flow Graph
  (CFG) and an **analysis object**, and perform analysis on the CFG.

There are also two subdirectories `apron` and `elina` that show how to
use Apron and Elina libraries. Currently, they are not compatible with
each other so that's why they must be built separately.
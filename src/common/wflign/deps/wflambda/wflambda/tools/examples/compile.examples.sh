#!/bin/bash

# Check library
if [ ! -f "../../build/libwflambda.a" ]; then
  echo "Library libwflambda.a not found. Please compile WFλ library from top folder first"
  exit 
fi

# Compile examples
g++ -O3 -I../.. -L../../build wfa_basic.cpp -o wfa_basic -lwflambda
g++ -O3 -I../.. -L../../build wfa_repeated.cpp -o wfa_repeated -lwflambda
g++ -O3 -I../.. -L../../build wfa_adapt.cpp -o wfa_adapt -lwflambda

#!/bin/bash

# Check library
if [ ! -f "../../build/libwfa.a" ]; then
  echo "Library libwfa.a not found. Please compile WFA library from top folder first"
  exit 
fi

# Compile examples
WFA_INCLUDES="../.."
WFA_LIB_PATH="../../build"
FLAGS="-O3"

gcc $FLAGS -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_basic.c -o wfa_basic -lwfa
gcc $FLAGS -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_repeated.c -o wfa_repeated -lwfa
gcc $FLAGS -I $WFA_INCLUDES -L $WFA_LIB_PATH wfa_adapt.c -o wfa_adapt -lwfa

WFA_CPP="../../bindings/cpp/WFAligner.cpp" 
g++ $FLAGS -I $WFA_INCLUDES -L $WFA_LIB_PATH $WFA_CPP wfa_bindings.cpp -o wfa_bindings -lwfa
g++ $FLAGS -I $WFA_INCLUDES -L $WFA_LIB_PATH $WFA_CPP wfa_lambda.cpp -o wfa_lambda -lwfa

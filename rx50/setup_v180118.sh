#!/bin/bash

LD_LIBRARY_PATH=/usr/local/cuda/lib64

SRC_DIR=$PWD/src
BIN_DIR=$PWD/bin

### rx50 software v180118
g++ -Wall $SRC_DIR/levelmonitor_v170215.cpp -o $BIN_DIR/levelmonitor -luhd -lboost_system -lboost_thread
g++ -Wall -c $SRC_DIR/sample_v180118.cpp -lpthread -luhd -lboost_system -lboost_thread
nvcc $PWD/sample_v180118.o $SRC_DIR/rx50_v170215.cu -o $BIN_DIR/rx50 -lboost_system -lcublas -lpthread -lcufft -lm -luhd -lgsl -lgslcblas
gcc -Wall $SRC_DIR/obs2ITU_v180118.c -o $BIN_DIR/obs2ITU -lm

#!/bin/bash

MX=`lscpu | grep '^CPU(s):' | awk '{print $2}'`
THREADS=`printf %.0f $(echo "${MX} * 0.93" | bc)`

sed -i "s/TASKSIZE :=.*;/TASKSIZE := ${THREADS};/" config.mpl

#!/bin/bash

# Usage : ./run-profiling
# This script runs all testcases for ms using gprof and perf

set -euo pipefail

# initialization

mkdir -p profiling/gprof 
mkdir -p profiling/perf

ROOT_DIR="../.."

# gprof

gcc -O2 -pg -o "$ROOT_DIR/ms-gprof" "$ROOT_DIR/ms.c" "$ROOT_DIR/streec.c" "$ROOT_DIR/rand1.c" -lm 

for i in {1..6}; do

    TESTCASE="params.case.0$i" 

    source $TESTCASE 

    "$ROOT_DIR/ms-gprof" $cmd > /dev/null

    gprof "$ROOT_DIR/ms-gprof" gmon.out > profiling/gprof/reporte-case-0$i.txt

done

rm gmon.out

#perf 

gcc -O2 -g -o "$ROOT_DIR/ms-perf" "$ROOT_DIR/ms.c" "$ROOT_DIR/streec.c" "$ROOT_DIR/rand1.c" -lm

for i in {1..6}; do

    TESTCASE="params.case.0$i" 

    source $TESTCASE 

    perf stat -r 3 -o "profiling/perf/stat-case-0$i.txt" -e cycles,instructions,branches,branch-misses,L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses "$ROOT_DIR/ms-perf" $cmd > /dev/null
    
    perf record -o "profiling/perf/perf-data-case-0$i.data" --call-graph dwarf,65528 "$ROOT_DIR/ms-perf" $cmd > /dev/null
        
    perf report -i "profiling/perf/perf-data-case-0$i.data" --stdio --children > "profiling/perf/perf-dwarf-case-0$i.txt"

done

rm "$ROOT_DIR/ms-gprof" 
rm "$ROOT_DIR/ms_perf"
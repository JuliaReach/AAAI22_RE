#!/bin/bash

make benchmarks

echo "Unicycle benchmark"
LD_LIBRARY_PATH=$HOME/programs/links/ ./B_Unicycle

echo "TORA benchmark"
LD_LIBRARY_PATH=$HOME/programs/links/ ./B_TORA

echo "ACC benchmark"
LD_LIBRARY_PATH=$HOME/programs/links/ ./B_ACC

echo "SinglePendulum benchmark"
LD_LIBRARY_PATH=$HOME/programs/links/ ./B_SinglePendulum

echo "DoublePendulum benchmark"
LD_LIBRARY_PATH=$HOME/programs/links/ ./B_DoublePendulum

echo "Airplane benchmark"
LD_LIBRARY_PATH=$HOME/programs/links/ ./B_Airplane

cd Results
./pandoc results_table -V geometry:margin=1in -V fontsize:1pt -s -o results_table.pdf

#!/bin/bash

sets='15 10 5 1'
for perc in $sets
do
echo SR"$perc"_DDT5_DATA
./BKGOnlyFit.sh SR"$perc"_DDT5_DATA 3 3
./GetSimpleLimit.sh SR"$perc"_DDT5_DATA_p3r3
./BKGOnlyFit.sh SR"$perc"_DDT5_DATA 3 3 --pseudo
./GetSimpleLimit.sh SR"$perc"_DDT5_DATA--pseudo_p3r3
done

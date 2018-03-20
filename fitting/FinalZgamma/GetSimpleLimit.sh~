#!/bin/bash

echo "------------------------"

mv OUTPUT_$1/card.txt .
mv OUTPUT_$1/base.root .
mv OUTPUT_$1/ralphabase.root .

cp templates/blanklimits.txt Limits_wC_$1.txt
sed -i 's/zqq50/SIGMASS/g' card.txt
masses='zqq10 zqq25 zqq50 zqq75 zqq100 zqq125'
for mass in $masses
do
sed -i 's/SIGMASS/'$mass'/g' card.txt
combine -M Asymptotic card.txt --rMin 0 --rMax 50
python GetLimitsForPoint.py Limits_wC_$1.txt
sed -i 's/'$mass'/SIGMASS/g' card.txt
done	

mv card.txt OUTPUT_$1/
mv base.root OUTPUT_$1/
mv ralphabase.root OUTPUT_$1/


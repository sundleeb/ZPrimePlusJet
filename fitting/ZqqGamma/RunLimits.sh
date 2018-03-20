#!/bin/bash

echo "------------------------"
cp templates/blanklimits.txt theselimits.txt
sed -i 's/zqq50/SIGMASS/g' cards_SRall.txt
#masses='zqq10 zqq15 zqq20 zqq25 zqq30 zqq35 zqq40 zqq45 zqq50 zqq55 zqq60 zqq65 zqq70 zqq75 zqq80 zqq85 zqq90 zqq95 zqq100 zqq105 zqq110 zqq115 zqq120 zqq125'
masses='zqq10 zqq25 zqq50 zqq75 zqq100 zqq125'
for mass in $masses
do
sed -i 's/SIGMASS/'$mass'/g' cards_SRall.txt
combine -M Asymptotic cards_SRall.txt --rMin -50 --rMax 50
python GetLimitsForPoint.py
sed -i 's/'$mass'/SIGMASS/g' cards_SRall.txt
done
#python AddToLimitPlot.py

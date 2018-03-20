
cp OUTPUT_$1/card.txt card_BIAS.txt
cp OUTPUT_$1/base.root .
cp OUTPUT_$1/ralphabase.root .

sed -i 's/zqq50/SIGMASS/g' card_BIAS.txt
masses='zqq10 zqq15 zqq20 zqq25 zqq30 zqq35 zqq50 zqq75 zqq100 zqq125'
for mass in $masses
do
echo $1$mass
sed -i 's/SIGMASS/'$mass'/g' card_BIAS.txt
python StatisticalTests.py $1$mass
sed -i 's/'$mass'/SIGMASS/g' card_BIAS.txt
done	


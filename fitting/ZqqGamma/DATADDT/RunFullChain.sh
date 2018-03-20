#!/bin/bash
echo "---- ---- ---- ----"
echo "treemaker"
#python MicroTreeProducer.py
#rm nonres.root
#hadd nonres.root gjets.root qcd.root
echo "DDT"
python DDTMaker.py
#python DDTClosure.py
echo "preproc SR"
python PrepareRalphabetData.py
echo "preproc SB"
python PrepareRalphabetData.py --T yes
cp  ZGammaFermiSRTemplates.root ..
cp  ZGammaFermiTTSBTemplates.root ..
ls



#!/bin/bash
echo "---- ---- ---- ----"
echo "treemaker"
python DataPreProc.py
rm nonres.root
hadd nonres.root gjets.root qcd.root
echo "DDT"
python DDTMaker.py
echo "preproc SR"
python PrepareRalphabetData.py
echo "preproc SB"
python PrepareRalphabetData.py --T yes
mv  ZGammaFermiSRTemplates.root ..
mv  ZGammaFermiTTSBTemplates.root ..
ls



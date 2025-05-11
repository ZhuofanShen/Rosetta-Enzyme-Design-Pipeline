for path in `ls HEM/*/cyclopropanation_styrene_EDA_distal/*_FastDesign/*rot*/*.sh`; do cd ${path%/*}; bash ${path##*/}; sleep 0.01; cd ../../../../..; done

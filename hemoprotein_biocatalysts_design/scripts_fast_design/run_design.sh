for path in `ls HEM/*/cyclopropanation_EDA_styrene/*/*/*rot*/design.sh`; do cd ${path%/*}; bash ${path##*/}; sleep 0.01; cd ../../../../../..; done

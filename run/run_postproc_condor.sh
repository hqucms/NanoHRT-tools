#!/bin/bash
#set -o verbose
workdir=`pwd`

echo `hostname`
echo "workdir: $workdir"
echo "args: $@"
ls -l

jobid=$1

source /cvmfs/cms.cern.ch/cmsset_default.sh
tar -xf CMSSW*.tar.gz --warning=no-timestamp

echo "printing pwd"
pwd

### --------------------------------###
#Keep track of release sandbox version
basedir=$PWD
rel=$(echo CMSSW_*)
arch=$(ls $rel/.SCRAM/|grep el) || echo "Failed to determine SL release!"
old_release_top=$(awk -F= '/RELEASETOP/ {print $2}' $rel/.SCRAM/el*/Environment) || echo "Failed to determine old releasetop!"
 
# Creating new release
# This is done so e.g CMSSW_BASE and other variables are not hardcoded to the sandbox setting paths
# which will not exist here
 
echo ">>> creating new release $rel"
mkdir tmp
cd tmp
export SCRAM_ARCH="$arch"
scramv1 project -f CMSSW $rel
new_release_top=$(awk -F= '/RELEASETOP/ {print $2}' $rel/.SCRAM/el*/Environment)
cd $rel
echo ">>> preparing sandbox release $rel"
 
for i in biglib bin cfipython config external include lib python src; do
    rm -rf "$i"
    mv "$basedir/$rel/$i" .
done
 
 
echo ">>> fixing python paths"
for f in $(find -iname __init__.py); do
    sed -i -e "s@$old_release_top@$new_release_top@" "$f"
done
 
eval $(scramv1 runtime -sh) || echo "The command 'cmsenv' failed!"
cd "$basedir"
echo "[$(date '+%F %T')] wrapper ready"
### --------------------------------###

ls -l

export MLAS_DYNAMIC_CPU_ARCH=99
export TMPDIR=`pwd`
python3 processor.py $jobid
status=$?

ls -l

if [ $status -ne 0 ]; then
    rm *.root
    ls -l
fi
#set +o verbose
exit $status

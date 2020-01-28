# NanoHRT-tools

### Set up CMSSW

```bash
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
cmsenv
```

### Set up official NanoAOD-tools

```bash
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
```

### Get customized NanoAOD tools for HeavyResTagging (NanoHRT-tools)

```bash
git clone https://github.com/hqucms/NanoHRT-tools.git PhysicsTools/NanoHRTTools
```

### Compile

```bash
scram b -j8
```

### Test

Instructions to run the nanoAOD postprocessor can be found at [nanoAOD-tools](https://github.com/cms-nanoAOD/nanoAOD-tools#nanoaod-tools). 

### Production

```bash
cd PhysicsTools/NanoHRTTools/run
```

 - To make trees for MC performance study:

```bash
python runPostProcessing.py -i /path/of/input -o /path/to/output --friend -I PhysicsTools.NanoHRTTools.producers.hrtMCTreeProducer hrtMCTree
```

 - To make trees for data/MC comparison and scale factor measurement:

```bash
python runHRTTrees.py -i /eos/uscms/store/user/lpcjme/noreplica/NanoHRT/path/to/input -o /path/to/output --channel 
[muon|photon|qcd] --year [2016|2017|2018] -n 20 --batch
```

  - the preselection for each channel is coded in `runHRTTrees.py`
  - add `--run-data` to make data trees
  - add `--run-syst` to make the systematic trees
  - the `--batch` option will submit jobs to condor automatically wihtout confirmation
     
More options of `runPostProcessing.py` or `runHRTTrees.py` (a wrapper of `runPostProcessing.py`) can be found with `python runPostProcessing.py -h` or `python runHRTTrees.py -h`, e.g.,

 - To resubmit failed jobs, run the same command but add `--resubmit`.

 - To add cross section weights and merge output trees according to the config file, run the same command but add `--post`.


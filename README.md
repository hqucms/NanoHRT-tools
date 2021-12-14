# NanoHRT-tools

### Set up CMSSW and offcial NanoAOD-tools

```bash
cmsrel CMSSW_11_1_0_pre5_PY3
cd CMSSW_11_1_0_pre5_PY3/src
cmsenv

git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
```

### Get customized NanoAOD tools for HeavyResTagging (NanoHRT-tools)

```bash
git clone https://github.com/hqucms/NanoHRT-tools.git PhysicsTools/NanoHRTTools -b dev/UL
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

##### Make trees for MC performance study:

```bash
python runPostProcessing.py [-i /path/of/input] -o /path/to/output -d datasets.yaml --friend 
-I PhysicsTools.NanoHRTTools.producers.hrtMCTreeProducer hrtMCTree -n 1
```

To merge the trees, run the same command but add `--post -w ''` (i.e., set `-w` to an empty string (`''`) -- we do not add the cross sections, but simply reweight signals to match the QCD spectrum afterwards).


##### Make trees for heavy flavour tagging (bb/cc) or top/W data/MC comparison and scale factor measurement:

```bash
python runHeavyFlavTrees.py -i /eos/uscms/store/user/lpcjme/noreplica/NanoHRT/path/to/input -o /path/to/output 
(--sample-dir custom_samples) --jet-type [ak8,ak15] --channel [photon|qcd|muon|inclusive] --year [2016|2017|2018] -n 10 
(--batch) (--run-data) (--run-syst)
(--run-tagger) (--run-mass-regression) (--sfbdt 0.5)
(--condor-extras '+AccountingGroup = "group_u_CMST3.all"')
```

Command line options:

  - the preselection for each channel is coded in `runHRTTrees.py`
  - add `--run-data` to make data trees
  - add `--run-syst` to make the systematic trees
  - can run data & MC for multiple years together w/ e.g., `--year 2016,2017,2018`. The `--run-data` option will be ignored in this case. Add also `--run-syst` to make the systematic trees.
  - use `--sample-dir` to specify the directory containing the sample lists. Currently we maintain two sets of sample lists: the default one is under [samples](run/samples) which is used for running over official NanoAOD datasets remotely, and the other one is [custom_samples](run/custom_samples) which is used for running over privately produced NanoAOD datasets locally. To run over the private produced samples, ones needs to add `--sample-dir custom_samples` to the command line.
  - the `--batch` option will submit jobs to condor automatically without confirmation
  - remove `-i` to run over remote files (e.g., official NanoAOD, or private NanoAOD published on DAS); consider adding `--prefetch` to copy files first before running
  - **[NEW]** add `--run-tagger` (`--run-mass-regression`) to run new ParticleNet tagger (mass regression) on-the-fly. Check `HeavyFlavBaseProducer.py` for the model configuration.
  - **[NEW]** use `--sfbdt` to change the sfBDT cut value. This affects only QCD and photon samples. By default, sfBDT > 0.5 is applied to QCD and photon samples.
  - **[NEW]** use `--condor-extras` to pass extra options to condor job description file.
     
More options of `runPostProcessing.py` or `runHRTTrees.py` (a wrapper of `runPostProcessing.py`) can be found with `python runPostProcessing.py -h` or `python runHRTTrees.py -h`, e.g.,

 - To resubmit failed jobs, run the same command but add `--resubmit`.

 - To add cross section weights and merge output trees according to the config file, run the same command but add `--post`. The cross section file to use can be set with the `-w` option. 


### Truth-matching criteria

For maximal flexibility, a number of truth-matching varibles are defined in [HeavyFlavBaseProducer](python/producers/HeavyFlavBaseProducer.py) for hadronically decaying top quarks and W, Z, Higgs bosons. For W/Z/H we define:

 - `fj_idx_dr_X`: deltaR of the fatjet to the nearest **hadronically decaying** X particle. **If found, this top quark `X` is then used to define all the following variables.** Default to 99 if no hadronically decaying X in the event.
 - `fj_idx_dr_X_daus`: max deltaR between the fatjet and the two quarks from X decay.
 - `fj_idx_X_pt`: pt of X
 - `fj_idx_X_decay`: max abs(pdgId) of the two quarks from X decay. For H/Z, this means 5: bb, 4: cc, <4: qq. For W, this means 4: cx, <4: qq. Default to 0 if no hadronically decaying X in the event.

Top quark is treated a bit differently:

 - `fj_idx_dr_T`: deltaR of the fatjet to the nearest **hadronically decaying** top quark. **If found, this top quark `T` is then used to define all the following variables.** Default to 99 if no hadronically decaying top in the event.
 - `fj_idx_dr_T_b`: deltaR between the fatjet and the b quark from the hadronic `T` decay.
 - `fj_idx_dr_T_Wq_(max|min)`: max|min deltaR between the fatjet and the two quarks from the W decay.
 - `fj_idx_T_Wq_(max|min)_pdgId`: pdgId (**w/o** taking the absolute value) of the corresponding two quarks from W decay.
 - `fj_idx_T_pt`: pt of `T`

#### Truth-matching criteria for top/W tagging scale factors


 - top-matched: all three quarks contained in the fatjet
   - `fj_1_dr_T_b<jetR && fj_1_dr_T_Wq_max<jetR`
 - W-matched: only the two W quarks contained, the b quark is outside the jet cone (if the W is from top quark decay)
   - `((fj_1_T_Wq_max_pdgId==0 && fj_1_dr_W_daus<jetR) || (fj_1_T_Wq_max_pdgId!=0 && fj_1_dr_T_b>=jetR && fj_1_dr_T_Wq_max<jetR))`
   - **[Note]** the first part is mainly intended for tW events where the top quark decays leptonically, and the W boson decays hadronically. This can be a sizeable contribution to the W-matched events and needs to be taken into account properly. The trick here makes use of the fact that `fj_1_T_Wq_max_pdgId` is non-zero only if there is a hadronic top in the event.
 - unmatched: defined as `(NOT top-matched) and (NOT W-matched)`, i.e.,
   - `!(fj_1_dr_T_b<jetR && fj_1_dr_T_Wq_max<jetR) && !((fj_1_T_Wq_max_pdgId==0 && fj_1_dr_W_daus<jetR) || (fj_1_T_Wq_max_pdgId!=0 && fj_1_dr_T_b>=jetR && fj_1_dr_T_Wq_max<jetR))`

[Extra] For selecting specifically W->cx decays from the W-matched jets:

 - W(cx)-matched:
   - `((fj_1_T_Wq_max_pdgId==0 && fj_1_dr_W_daus<jetR && fj_1_W_decay==4) || (fj_1_T_Wq_max_pdgId!=0 && fj_1_dr_T_b>=jetR && fj_1_dr_T_Wq_max<jetR && (abs(fj_1_T_Wq_max_pdgId)==4 || abs(fj_1_T_Wq_min_pdgId)==4)))`

### Checklist when updating to new data-taking years / production campaigns

- [ ] triggers
- [ ] lumi values
- [ ] golden JSON
- [ ] PU rewgt
- [ ] lepton ID/ISO
- [ ] b-tag WP
- [ ] JEC/JER
- [ ] MET filters
- [ ] MET recipes (if any)
- [ ] samples (check also those in PRODUCTION status)

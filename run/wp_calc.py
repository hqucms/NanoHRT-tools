import ROOT
import uproot
import matplotlib.pyplot as plt
import numpy as np
import sys

era = str(sys.argv[1])

qcd_file = uproot.open("/eos/user/h/hlarson/NanoHRT/_ak8_qcd_" + era + "/mc/mg-qcd_tree.root")
tree = qcd_file["Events"]

print("we have a tree")

branches = tree.arrays()
print("we have branches")
xsecW = branches["xsecWeight"]
total = np.sum(xsecW)
print("we have a total : ", total)

qcd_taggers = []
for item in tree.keys():
    if item[0:16] == "fj_1_ParticleNet" and item[-3:] == "QCD":
        print("qcd_tagger added")
        qcd_taggers.append(item)
fraction = []


with ROOT.TFile("qcd_wp_histos.root", "recreate") as outfile:
    for tagger in qcd_taggers:
        tag_type = tagger.split("_")[-1]
        histo = ROOT.TH1D(name=f"qcd_{tag_type}", title=tag_type, nbinsx=600, xlow=0.4, xup=1.0)
        print("histo initialized")

        for counter in np.arange(0.3, 1, 0.001):
            passing = 0
            for idx, item in enumerate(branches[tagger]):
                if item > counter:
                    passing += xsecW[idx]
            fraction = passing/total

            histo.Fill(counter, fraction)

        print("we have a fraction! heres the last one: ", fraction)
        histo.Draw()

        outfile.WriteObject(histo, histo.GetName())

print("Done with wp :D")
    

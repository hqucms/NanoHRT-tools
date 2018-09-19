import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

def get_corrected_sdmass(jet, subjets):
    if not jet or not subjets or len(subjets)==0:
        return 0
    msd_uncorr = sum([sj.p4() * (1 - sj.rawFactor) for sj in subjets], ROOT.TLorentzVector()).M()

    pt = jet.pt
    eta = jet.eta
    gencorr = 1.006261 + ((-1.061605) * pow(pt*0.079990,-1.204538))
    recocorr = 1
    if abs(eta) <= 1.3:
        recocorr = 1.093020 + (-0.000150068) * pt + (3.44866e-07) * pow(pt, 2) + (-2.68100e-10) * pow(pt, 3) + (8.67440e-14) * pow(pt, 4) + (-1.00114e-17) * pow(pt, 5)
    else:
        recocorr = 1.272115 + (-0.000571640) * pt + (8.37289e-07) * pow(pt, 2) + (-5.20433e-10) * pow(pt, 3) + (1.45375e-13) * pow(pt, 4) + (-1.50389e-17) * pow(pt, 5)

    return msd_uncorr * gencorr * recocorr


def get_sdmass_fromsubjets(jet, subjets):
    if not jet or not subjets or len(subjets) == 0:
        return 0
    msd = sum([sj.p4() for sj in subjets], ROOT.TLorentzVector()).M()
    return msd


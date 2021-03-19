import math
import logging
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


def clip(value, lower, upper):
    return lower if value < lower else upper if value > upper else value


def deltaPhi(phi1, phi2):
    try:
        dphi = phi1 - phi2
    except TypeError:
        dphi = phi1.phi - phi2.phi
    while dphi > math.pi:
        dphi -= 2 * math.pi
    while dphi < -math.pi:
        dphi += 2 * math.pi
    return dphi


def deltaR2(eta1, phi1, eta2=None, phi2=None):
    if eta2 is None:
        a, b = eta1, phi1
        return deltaR2(a.eta, a.phi, b.eta, b.phi)
    else:
        deta = eta1 - eta2
        dphi = deltaPhi(phi1, phi2)
        return deta * deta + dphi * dphi


def deltaR(eta1, phi1, eta2=None, phi2=None):
    return math.sqrt(deltaR2(eta1, phi1, eta2, phi2))


def closest(obj, collection, presel=lambda x, y: True):
    ret = None
    dr2Min = 1e6
    for x in collection:
        if not presel(obj, x):
            continue
        dr2 = deltaR2(obj, x)
        if dr2 < dr2Min:
            ret = x
            dr2Min = dr2
    return (ret, math.sqrt(dr2Min))


def polarP4(obj=None, pt='pt', eta='eta', phi='phi', mass='mass'):
    if obj is None:
        return ROOT.Math.PtEtaPhiMVector()
    pt_val = getattr(obj, pt) if pt else 0
    eta_val = getattr(obj, eta) if eta else 0
    phi_val = getattr(obj, phi) if phi else 0
    mass_val = getattr(obj, mass) if mass else 0
    return ROOT.Math.PtEtaPhiMVector(pt_val, eta_val, phi_val, mass_val)


def p4(obj=None, pt='pt', eta='eta', phi='phi', mass='mass'):
    v = polarP4(obj, pt, eta, phi, mass)
    return ROOT.Math.XYZTVector(v.px(), v.py(), v.pz(), v.energy())


def sumP4(*args):
    p4s = [polarP4(x) for x in args]
    return sum(p4s, ROOT.Math.PtEtaPhiMVector())


def p4_str(p):
    return '(pt=%s, eta=%s, phi=%s, mass=%s)' % (p.pt(), p.eta(), p.phi(), p.mass())


def get_subjets(jet, subjetCollection, idxNames=('subJetIdx1', 'subJetIdx2')):
    subjets = []
    for idxname in idxNames:
        idx = getattr(jet, idxname)
        if idx >= 0:
            subjets.append(subjetCollection[idx])
    subjets = sorted(subjets, key=lambda x: x.pt, reverse=True)  # sort by pt
    return subjets


def corrected_svmass(sv):
    pproj = polarP4(sv).P() * math.sin(sv.pAngle)
    return math.sqrt(sv.mass * sv.mass + pproj * pproj) + pproj


def transverseMass(obj, met):
    cos_dphi = math.cos(deltaPhi(obj, met))
    return math.sqrt(2 * obj.pt * met.pt * (1 - cos_dphi))


def minValue(collection, fallback=99):
    if len(collection) == 0:
        return fallback
    else:
        return min(collection)


def maxValue(collection, fallback=0):
    if len(collection) == 0:
        return fallback
    else:
        return max(collection)


def configLogger(name, loglevel=logging.INFO, filename=None):
    # define a Handler which writes INFO messages or higher to the sys.stderr
    logger = logging.getLogger(name)
    logger.setLevel(loglevel)
    console = logging.StreamHandler()
    console.setLevel(loglevel)
    console.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s'))
    logger.addHandler(console)
    if filename:
        logfile = logging.FileHandler(filename)
        logfile.setLevel(loglevel)
        logfile.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s'))
        logger.addHandler(logfile)

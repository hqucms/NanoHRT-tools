from .HeavyFlavPhotonSampleProducer import PhotonSampleProducer
from .HeavyFlavQCDSampleProducer import QCDSampleProducer
from .HeavyFlavMuonSampleProducer import MuonSampleProducer
from .HeavyFlavInclusiveSampleProducer import InclusiveSampleProducer


def heavyFlavSFTreeFromConfig():
    import yaml
    with open('heavyFlavSFTree_cfg.json') as f:
        cfg = yaml.safe_load(f)
        channel = cfg['channel']
        del cfg['channel']
    if channel == 'photon':
        return PhotonSampleProducer(**cfg)
    elif channel == 'qcd':
        return QCDSampleProducer(**cfg)
    elif channel == 'muon':
        return MuonSampleProducer(**cfg)
    elif channel == 'inclusive':
        return InclusiveSampleProducer(**cfg)
    else:
        return RuntimeError('Unsupported channel %s' % channel)

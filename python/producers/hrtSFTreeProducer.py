from PhysicsTools.NanoHRTTools.producers.MuonSampleProducer import MuonSampleProducer
from PhysicsTools.NanoHRTTools.producers.PhotonSampleProducer import PhotonSampleProducer
from PhysicsTools.NanoHRTTools.producers.QCDSampleProducer import QCDSampleProducer


def hrtSFTreeFromConfig():
    import yaml
    with open('hrtSFTree_cfg.json') as f:
        cfg = yaml.safe_load(f)
        channel = cfg['channel']
        del cfg['channel']
    if channel == 'muon':
        return MuonSampleProducer(**cfg)
    elif channel == 'photon':
        return PhotonSampleProducer(**cfg)
    elif channel == 'qcd':
        return QCDSampleProducer(**cfg)
    else:
        return RuntimeError('Unsupported channel %s' % channel)


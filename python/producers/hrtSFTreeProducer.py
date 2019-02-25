from PhysicsTools.NanoHRTTools.producers.MuonSampleProducer import MuonSampleProducer
from PhysicsTools.NanoHRTTools.producers.PhotonSampleProducer import PhotonSampleProducer
from PhysicsTools.NanoHRTTools.producers.QCDSampleProducer import QCDSampleProducer


def hrtSFTreeFromConfig():
    import yaml
    with open('hrtSFTree_cfg.json') as f:
        cfg = yaml.safe_load(f)
    if cfg['channel'] == 'muon':
        return MuonSampleProducer(**cfg)
    elif cfg['channel'] == 'photon':
        return PhotonSampleProducer(**cfg)
    elif cfg['channel'] == 'qcd':
        return QCDSampleProducer(**cfg)
    else:
        return RuntimeError('Unsupported channel %s' % cfg['channel'])


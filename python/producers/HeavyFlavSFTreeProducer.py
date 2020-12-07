from PhysicsTools.NanoHRTTools.producers.HeavyFlavPhotonSampleProducer import PhotonSampleProducer
from PhysicsTools.NanoHRTTools.producers.HeavyFlavQCDSampleProducer import QCDSampleProducer
from PhysicsTools.NanoHRTTools.producers.HeavyFlavSignalSampleProducer import SignalSampleProducer
<<<<<<< HEAD
=======
from PhysicsTools.NanoHRTTools.producers.HeavyFlavInclusiveSampleProducer import InclusiveSampleProducer

>>>>>>> upstream/master

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
    elif channel == 'signal':
        return SignalSampleProducer(**cfg)
    elif channel == 'inclusive':
        return InclusiveSampleProducer(**cfg)
    else:
        return RuntimeError('Unsupported channel %s' % channel)


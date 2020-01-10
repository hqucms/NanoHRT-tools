from PhysicsTools.NanoHRTTools.producers.HFPhotonSampleProducer import HFPhotonSFTreeProducer
from PhysicsTools.NanoHRTTools.producers.HFQCDSampleProducer import HFQCDSFTreeProducer


def hrtSFTreeFromConfig():
    import yaml
    with open('hrtSFTree_cfg.json') as f:
        cfg = yaml.safe_load(f)
        channel = cfg['channel']        
        del cfg['channel']   
    if (channel == 'photon'):
	return HFPhotonSFTreeProducer(**cfg)
    elif(channel == 'qcd'):
        return HFQCDSFTreeProducer(**cfg)

    else:
        return RuntimeError('Unsupported channel %s' % channel)

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from .HeavyFlavBaseProducer import HeavyFlavBaseProducer


class QCDSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(QCDSampleProducer, self).__init__(channel='qcd', **kwargs)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(QCDSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        # trigger variables
        self.out.branch("passHTTrig", "O")

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.selectLeptons(event)
        self.correctJetsAndMET(event)

        if len(event.fatjets) < 2:
            return False
        probe_jets = event.fatjets[:2]

        if self._opts['sfbdt_threshold'] > 0:
            self.selectSV(event)
            if event.secondary_vertices < 2:
                return False
            self.matchSVToFatJets(event, probe_jets)
            # require: 2 subjets, 50 < sdmass < 200, sfBDT > 0.5 (and >=1 SV matched to each subjet -- included in the sfBDT requirement)
            probe_jets = [fj for fj in probe_jets
                          if len(fj.subjets) == 2 and fj.msoftdrop > 50 and fj.msoftdrop < 200
                          and fj.sfBDT > self._opts['sfbdt_threshold']]
            if len(probe_jets) == 0:
                return False

        self.loadGenHistory(event, probe_jets)
        self.evalTagger(event, probe_jets)
        self.evalMassRegression(event, probe_jets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probe_jets)

        if self.year == 2016:
            self.out.fillBranch("passHTTrig", event.HLT_PFHT900)
        else:
            self.out.fillBranch("passHTTrig", event.HLT_PFHT1050)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def QCDTree_2016(): return QCDSampleProducer(year=2016)
def QCDTree_2017(): return QCDSampleProducer(year=2017)
def QCDTree_2018(): return QCDSampleProducer(year=2018)

import numpy as np
import awkward
import uproot
from uproot_methods import TLorentzVectorArray
from collections import Counter


class ParticleNetTagInfoMaker(object):

    def __init__(self, fatjet_branch='FatJet', pfcand_branch='PFCands', sv_branch='SV', jetR=0.8, pfcand_ptcut=0):
        self.fatjet_branch = fatjet_branch
        self.pfcand_branch = pfcand_branch
        self.sv_branch = sv_branch
        self.idx_branch = '{jet}To{cand}_candIdx'.format(jet=fatjet_branch, cand=pfcand_branch)
        self.jet_r2 = jetR * jetR
        self.pfcand_ptcut = pfcand_ptcut

    def _finalize_data(self, data):
        for k in data:
            data[k] = data[k].astype('float32')

    def _make_pfcands(self, table):
        data = {}
        if self.idx_branch in table:
            jet_cand_counts = table[self.fatjet_branch + '_nPFCand']
            jet_cand_idxs = table[self.idx_branch]
        else:
            cand_parents = table[self.pfcand_branch + '_jetIdx']
            c = Counter((cand_parents.offsets[:-1] + cand_parents).content)
            jet_cand_counts = awkward.JaggedArray.fromcounts(self.jetp4.counts, [c[k] for k in sorted(c.keys())])

        ptcut = None

        def pf(var_name):
            branch_name = self.pfcand_branch + '_' + var_name
            if var_name[:4] == 'btag' and self.idx_branch in table:
                branch_name = branch_name + '_' + self.fatjet_branch
            cand_arr = table[branch_name]
            if self.idx_branch in table:
                cand_arr = cand_arr[jet_cand_idxs]
            out = jet_cand_counts.copy(
                content=awkward.JaggedArray.fromcounts(jet_cand_counts.content, cand_arr.content))
            if ptcut is None:
                return out
            else:
                return out[ptcut]

        if self.pfcand_ptcut > 0:
            ptcut = pf('pt') > self.pfcand_ptcut

        data['pfcand_VTX_ass'] = pf('pvAssocQuality')
        data['pfcand_lostInnerHits'] = pf('lostInnerHits')
        data['pfcand_quality'] = pf('trkQuality')

        pdgId = pf('pdgId')
        charge = pf('charge')
        data['pfcand_isEl'] = np.abs(pdgId) == 11
        data['pfcand_isMu'] = np.abs(pdgId) == 13
        data['pfcand_isChargedHad'] = np.abs(pdgId) == 211
        data['pfcand_isGamma'] = np.abs(pdgId) == 22
        data['pfcand_isNeutralHad'] = np.abs(pdgId) == 130
        data['pfcand_charge'] = charge

        dz = pf('dz')
        dxy = pf('d0')
#         ### FIXME ###
#         is_neutral = np.logical_or(charge.content == 0, np.logical_and(dz.content == -1.0, dxy.content == -1.0))
#         dz.content[is_neutral] = 0
#         dxy.content[is_neutral] = 0
#         ### FIXME ###
        data['pfcand_dz'] = dz
        data['pfcand_dxy'] = dxy
        data['pfcand_dzsig'] = dz / pf('dzErr')
        data['pfcand_dxysig'] = dxy / pf('d0Err')

        candp4 = TLorentzVectorArray.from_ptetaphim(pf('pt'), pf('eta'), pf('phi'), pf('mass'))
        data['pfcand_mask'] = charge.ones_like()
        data['pfcand_phirel'] = candp4.delta_phi(self.jetp4)
        data['pfcand_etarel'] = self.eta_sign * (candp4.eta - self.jetp4.eta)
        data['pfcand_abseta'] = np.abs(candp4.eta)

        data['pfcand_pt_log_nopuppi'] = np.log(candp4.pt)
        data['pfcand_e_log_nopuppi'] = np.log(candp4.energy)

        chi2 = pf('trkChi2')
        chi2.content[chi2.content == -1] = 999
        data['pfcand_normchi2'] = np.floor(chi2)
        data['pfcand_btagEtaRel'] = pf('btagEtaRel')
        data['pfcand_btagPtRatio'] = pf('btagPtRatio')
        data['pfcand_btagPParRatio'] = pf('btagPParRatio')
        data['pfcand_btagSip3dVal'] = pf('btagSip3dVal')
        data['pfcand_btagSip3dSig'] = pf('btagSip3dSig')
        data['pfcand_btagJetDistVal'] = pf('btagJetDistVal')

        self._finalize_data(data)
        self.data.update(data)

    def _make_sv(self, table):
        data = {}
        all_svp4 = TLorentzVectorArray.from_ptetaphim(
            table[self.sv_branch + '_pt'],
            table[self.sv_branch + '_eta'],
            table[self.sv_branch + '_phi'],
            table[self.sv_branch + '_mass'],
        )

        jet_cross_sv = self.jetp4.cross(all_svp4, nested=True)
        match = jet_cross_sv.i0.delta_r2(jet_cross_sv.i1) < self.jet_r2

        def sv(var_name):
            sv_arr = table[self.sv_branch + '_' + var_name]
            return self.jetp4.eta.cross(sv_arr, nested=True).unzip()[1][match]

        svp4 = TLorentzVectorArray.from_ptetaphim(sv('pt'), sv('eta'), sv('phi'), sv('mass'))
        data['sv_phirel'] = svp4.delta_phi(self.jetp4)
        data['sv_etarel'] = self.eta_sign * (svp4.eta - self.jetp4.eta)
        data['sv_abseta'] = np.abs(svp4.eta)
        data['sv_mass'] = svp4.mass
        data['sv_pt_log'] = np.log(svp4.pt)
        data['sv_mask'] = data['sv_pt_log'].ones_like()

        data['sv_ntracks'] = sv('ntracks')
        data['sv_normchi2'] = sv('chi2')
        data['sv_dxy'] = sv('dxy')
        data['sv_dxysig'] = sv('dxySig')
        data['sv_d3d'] = sv('dlen')
        data['sv_d3dsig'] = sv('dlenSig')
        data['sv_costhetasvpv'] = -np.cos(sv('pAngle'))

        dxysig = sv('dxySig')
        dxysig.content[~np.isfinite(dxysig.content)] = 0
        pos = dxysig.argsort()
        for k in data:
            data[k] = data[k][pos]

        self._finalize_data(data)
        self.data.update(data)

    def convert(self, table):
        self.data = {}
        self.jetp4 = TLorentzVectorArray.from_ptetaphim(
            table[self.fatjet_branch + '_pt'],
            table[self.fatjet_branch + '_eta'],
            table[self.fatjet_branch + '_phi'],
            table[self.fatjet_branch + '_mass'],
        )
        self.eta_sign = self.jetp4.eta.ones_like()
        self.eta_sign[self.jetp4.eta <= 0] = -1
        self._make_pfcands(table)
        self._make_sv(table)
        self.data['_jetp4'] = self.jetp4
        return self.data

    def init_file(self, inputFile, fetch_step=1000):
        self._uproot_basketcache = uproot.cache.ThreadSafeArrayCache('500MB')
        self._uproot_keycache = uproot.cache.ThreadSafeArrayCache('10MB')
        self._uproot_tree = uproot.open(inputFile.GetName())['Events']
        self._uproot_fetch_step = fetch_step
        self._uproot_start = 0
        self._uproot_stop = 0
        self._taginfo = None

    def load(self, event_idx):
        if event_idx >= self._uproot_stop:
            # needs to fetch next batch
            self._uproot_start = event_idx
            self._uproot_stop = self._uproot_start + self._uproot_fetch_step
            table = self._uproot_tree.arrays(
                ['AK15PuppiToPFCands_candIdx', 'AK15Puppi_nPFCand', 'AK15Puppi_pt',
                 'AK15Puppi_eta', 'AK15Puppi_phi', 'AK15Puppi_mass', 'PFCands*', 'SV*'],
                namedecode='utf-8', entrystart=self._uproot_start, entrystop=self._uproot_stop,
                basketcache=self._uproot_basketcache, keycache=self._uproot_keycache,
            )
            self._taginfo = self.convert(table)
        return self._taginfo


if __name__ == '__main__':
    import uproot
    import argparse
    parser = argparse.ArgumentParser('TEST')
    parser.add_argument('-i', '--input')
    args = parser.parse_args()

    p = ParticleNetTagInfoMaker()
    table = uproot.open(args.input)['Events'].arrays(['FatJet*', 'PFCands*', 'SV*'], namedecode='utf-8')
    taginfo = p.convert(table)
    for k in taginfo:
        print(k, taginfo[k])

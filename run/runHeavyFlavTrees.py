#!/usr/bin/env python
from __future__ import print_function

import os
import copy

from runPostProcessing import get_arg_parser, run, tar_cmssw
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')

hrt_cfgname = 'heavyFlavSFTree_cfg.json'
default_config = {'sfbdt_threshold': -99,
                  'run_tagger': False, 'tagger_versions': ['V02b', 'V02c', 'V02d'],
                  'run_mass_regression': False, 'mass_regression_versions': ['V01a', 'V01b', 'V01c'],
                  'jec': False, 'jes': None, 'jes_source': '', 'jes_uncertainty_file_prefix': '',
                  'jer': 'nominal', 'jmr': None, 'met_unclustered': None, 'smearMET': False, 'applyHEMUnc': False}

cut_dict_ak8 = {
    'photon': 'Sum$(Photon_pt>200 && Photon_cutBased>=2 && Photon_electronVeto)>0 && nFatJet>0',
    'qcd': 'Sum$((Jet_pt>25 && abs(Jet_eta)<2.4 && (Jet_jetId & 2)) * Jet_pt)>200 && nFatJet>0',
    'muon': 'Sum$(Muon_pt>55 && abs(Muon_eta)<2.4 && Muon_tightId && Muon_miniPFRelIso_all<0.10)>0 && nFatJet>0',
    'diboson': '(Sum$(Electron_pt>20 && abs(Electron_eta)<2.5 && abs(Electron_dxy)<0.05 && abs(Electron_dz)<0.2 && Electron_mvaFall17V2noIso_WP90 && Electron_miniPFRelIso_all<0.4) >= 2 ||'
               ' Sum$(Muon_pt>20 && abs(Muon_eta)<2.4 && abs(Muon_dxy)<0.05 && abs(Muon_dz)<0.2 && Muon_looseId && Muon_miniPFRelIso_all<0.4) >= 2) && nFatJet>0',
    'inclusive': 'Sum$((Jet_pt>25 && abs(Jet_eta)<2.4 && (Jet_jetId & 2)) * Jet_pt)>300 && Sum$(FatJet_subJetIdx1>=0 && FatJet_subJetIdx2>=0 && FatJet_msoftdrop>10)>0',
}
cut_dict_ak15 = {
    'photon': 'Sum$(Photon_pt>200 && Photon_cutBased>=2 && Photon_electronVeto)>0 && nAK15Puppi>0',
    'qcd': 'Sum$((Jet_pt>25 && abs(Jet_eta)<2.4 && (Jet_jetId & 2)) * Jet_pt)>200 && nAK15Puppi>0',
    'muon': 'Sum$(Muon_pt>55 && abs(Muon_eta)<2.4 && Muon_tightId && Muon_miniPFRelIso_all<0.10)>0 && nAK15Puppi>0',
    'diboson': '(Sum$(Electron_pt>20 && abs(Electron_eta)<2.5 && abs(Electron_dxy)<0.05 && abs(Electron_dz)<0.2 && Electron_mvaFall17V2noIso_WP90 && Electron_miniPFRelIso_all<0.4) >= 2 ||'
               ' Sum$(Muon_pt>20 && abs(Muon_eta)<2.4 && abs(Muon_dxy)<0.05 && abs(Muon_dz)<0.2 && Muon_looseId && Muon_miniPFRelIso_all<0.4) >= 2) && nAK15Puppi>0',
    'inclusive': 'Sum$((Jet_pt>25 && abs(Jet_eta)<2.4 && (Jet_jetId & 2)) * Jet_pt)>300 && Sum$(AK15Puppi_subJetIdx1>=0 && AK15Puppi_subJetIdx2>=0 && AK15Puppi_msoftdrop>10)>0',
}

golden_json = {
    2015: 'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',
    2016: 'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',
    2017: 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',
    2018: 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt',
}


def _process(args):
    default_config['jetType'] = args.jet_type
    if args.run_tagger:
        default_config['run_tagger'] = True
        if args.jet_type == 'ak8':
            raise NotImplementedError('No training for ak8')
        logging.info('Will run tagger version(s): %s' % ','.join(default_config['tagger_versions']))
    if args.run_mass_regression:
        default_config['run_mass_regression'] = True
        if args.jet_type == 'ak8':
            default_config['mass_regression_versions'] = ['ak8V01a', 'ak8V01b', 'ak8V01c']
        logging.info('Will run mass regression version(s): %s' % ','.join(default_config['mass_regression_versions']))

    year = int(args.year)
    channel = args.channel
    default_config['year'] = year
    default_config['channel'] = channel
    if channel in ('qcd', 'photon'):
        default_config['sfbdt_threshold'] = args.sfbdt

    if year in (2017, 2018):
        args.weight_file = 'samples/xsec_2017.conf'

    basename = os.path.basename(args.outputdir) + '_' + args.jet_type + '_' + channel + '_' + str(year)
    args.outputdir = os.path.join(os.path.dirname(args.outputdir), basename, 'data' if args.run_data else 'mc')
    args.jobdir = os.path.join('jobs_%s' % basename, 'data' if args.run_data else 'mc')

    if args.run_data:
        args.datasets = '%s/%s_%d_DATA.yaml' % (args.sample_dir, channel, year)
        args.extra_transfer = os.path.expandvars(
            '$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/JSON/%s' % golden_json[year])
        args.json = golden_json[year]
    else:
        args.datasets = '%s/%s_%d_MC.yaml' % (args.sample_dir, channel, year)

    if args.jet_type == 'ak15':
        args.cut = cut_dict_ak15[channel]
    else:
        args.cut = cut_dict_ak8[channel]

    args.imports = [('PhysicsTools.NanoHRTTools.producers.HeavyFlavSFTreeProducer', 'heavyFlavSFTreeFromConfig')]
    if not args.run_data:
        args.imports.extend([('PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer',
                              'puWeight_UL2016' if year == 2015 else 'puWeight_UL%d' % year),
                             ('PhysicsTools.NanoHRTTools.producers.topPtWeightProducer', 'topPtWeight')])

    # data, or just nominal MC
    if args.run_data or not args.run_syst:
        cfg = copy.deepcopy(default_config)
        if args.run_data:
            cfg['jes'] = None
            cfg['jer'] = None
            cfg['jmr'] = None
            cfg['met_unclustered'] = None
        run(args, configs={hrt_cfgname: cfg})
        return

    # MC for syst.
    if args.run_syst and not args.run_data:

        # nominal w/ PDF/Scale wegihts
        logging.info('Start making nominal trees with PDF/scale weights...')
        syst_name = 'LHEWeight'
        opts = copy.deepcopy(args)
        cfg = copy.deepcopy(default_config)
        opts.outputdir = os.path.join(os.path.dirname(opts.outputdir), syst_name)
        opts.jobdir = os.path.join(os.path.dirname(opts.jobdir), syst_name)
        opts.branchsel_out = 'keep_and_drop_output_LHEweights.txt'
        run(opts, configs={hrt_cfgname: cfg})

        # JES up/down
        for variation in ['up', 'down']:
            syst_name = 'jes_%s' % variation
            logging.info('Start making %s trees...' % syst_name)
            opts = copy.deepcopy(args)
            cfg = copy.deepcopy(default_config)
            cfg['jes'] = variation
            opts.outputdir = os.path.join(os.path.dirname(opts.outputdir), syst_name)
            opts.jobdir = os.path.join(os.path.dirname(opts.jobdir), syst_name)
            run(opts, configs={hrt_cfgname: cfg})

        # JER up/down
        for variation in ['up', 'down']:
            syst_name = 'jer_%s' % variation
            logging.info('Start making %s trees...' % syst_name)
            opts = copy.deepcopy(args)
            cfg = copy.deepcopy(default_config)
            cfg['jer'] = variation
            opts.outputdir = os.path.join(os.path.dirname(opts.outputdir), syst_name)
            opts.jobdir = os.path.join(os.path.dirname(opts.jobdir), syst_name)
            run(opts, configs={hrt_cfgname: cfg})

        # MET unclustEn up/down
        if args.channel == 'muon':
            for variation in ['up', 'down']:
                syst_name = 'met_%s' % variation
                logging.info('Start making %s trees...' % syst_name)
                opts = copy.deepcopy(args)
                cfg = copy.deepcopy(default_config)
                cfg['met_unclustered'] = variation
                opts.outputdir = os.path.join(os.path.dirname(opts.outputdir), syst_name)
                opts.jobdir = os.path.join(os.path.dirname(opts.jobdir), syst_name)
                run(opts, configs={hrt_cfgname: cfg})


def main():
    parser = get_arg_parser()

    parser.add_argument('--jet-type',
                        choices=['ak8', 'ak15'],
                        required=True,
                        help='Jet type: ak8, ak15'
                        )

    parser.add_argument('--channel',
                        type=str,
                        required=True,
                        help='Channel: photon, qcd, muon, diboson, signal, inclusive, or comma separated list e.g., `qcd,photon`'
                        )

    parser.add_argument('--sfbdt',
                        type=float, default=0.5,
                        help='sfBDT cut, applies only to `qcd` and `photon` channels. Default: %(default)s'
                        )

    parser.add_argument('--run-syst',
                        action='store_true', default=False,
                        help='Run all the systematic trees. Default: %(default)s'
                        )

    parser.add_argument('--run-data',
                        action='store_true', default=False,
                        help='Run over data. Default: %(default)s'
                        )

    parser.add_argument('--year',
                        type=str,
                        required=True,
                        help='Year: 2015 (2016 preVFP), 2016 (2016 postVFP), 2017, 2018, or comma separated list e.g., `2016,2017,2018`'
                        )

    parser.add_argument('--sample-dir',
                        type=str,
                        default='samples',
                        help='Directory of the sample list files. Default: %(default)s'
                        )

    parser.add_argument('--run-tagger',
                        action='store_true', default=False,
                        help='Run tagger. Default: %(default)s'
                        )

    parser.add_argument('--run-mass-regression',
                        action='store_true', default=False,
                        help='Run mass regression. Default: %(default)s'
                        )

    args = parser.parse_args()

    if not (args.post or args.add_weight or args.merge):
        tar_cmssw(args.tarball_suffix)

    years = args.year.split(',')
    channels = args.channel.split(',')
    categories = ['data' if args.run_data else 'mc']

    for year in years:
        for chn in channels:
            for cat in categories:
                opts = copy.deepcopy(args)
                if cat == 'data':
                    opts.run_data = True
                    opts.nfiles_per_job *= 2
                if opts.inputdir:
                    opts.inputdir = opts.inputdir.rstrip('/').replace('_YEAR_', year)
                    assert(year in opts.inputdir)
                    if opts.inputdir.rsplit('/', 1)[1] not in ['data', 'mc']:
                        opts.inputdir = os.path.join(opts.inputdir, cat)
                    assert(opts.inputdir.endswith(cat))
                opts.year = year
                opts.channel = chn
                logging.info('inputdir=%s, year=%s, channel=%s, cat=%s, syst=%s', opts.inputdir, opts.year,
                             opts.channel, 'data' if opts.run_data else 'mc', 'syst' if opts.run_syst else 'none')
                _process(opts)


if __name__ == '__main__':
    main()

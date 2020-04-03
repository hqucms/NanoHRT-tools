#!/usr/bin/env python
from __future__ import print_function

import os
import copy
from runPostProcessing import get_arg_parser, run, tar_cmssw
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')

hrt_cfgname = 'heavyFlavSFTree_cfg.json'
default_config = {'data': False, 'jetType': None, 'channel': None, 'year': None, 'jec': False, 'jes': None, 'jes_source': '', 'jer': 'nominal', 'jmr': None, 'met_unclustered': None}
cut_dict_ak8 = {
    'photon': 'Sum$(Photon_pt>200 && (Photon_cutBasedBitmap & 2) && Photon_electronVeto)>0 && Sum$(FatJet_subJetIdx1>=0 && FatJet_subJetIdx2>=0 && FatJet_msoftdrop>10)>0',
    'qcd': 'Sum$((Jet_pt>25 && abs(Jet_eta)<2.4 && (Jet_jetId & 2)) * Jet_pt)>800 && Sum$(FatJet_subJetIdx1>=0 && FatJet_subJetIdx2>=0 && FatJet_msoftdrop>10)>0',
    'signal': 'Sum$(FatJet_subJetIdx1>=0 && FatJet_subJetIdx2>=0 && FatJet_msoftdrop>10)>0',
    }
cut_dict_ak15 = {
    'photon': 'Sum$(Photon_pt>200 && (Photon_cutBasedBitmap & 2) && Photon_electronVeto)>0 && Sum$(AK15Puppi_subJetIdx1>=0 && AK15Puppi_subJetIdx2>=0 && AK15Puppi_msoftdrop>10)>0',
    'qcd': 'Sum$((Jet_pt>25 && abs(Jet_eta)<2.4 && (Jet_jetId & 2)) * Jet_pt)>800 && Sum$(AK15Puppi_subJetIdx1>=0 && AK15Puppi_subJetIdx2>=0 && AK15Puppi_msoftdrop>10)>0',
    'signal': 'Sum$(AK15Puppi_subJetIdx1>=0 && AK15Puppi_subJetIdx2>=0 && AK15Puppi_msoftdrop>10)>0',
    }

golden_json = {
    2016: 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt',
    2017: 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt',
    2018: 'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt',
    }


def _process(args):
    default_config['jetType'] = args.jet_type

    channel = args.channel
    default_config['channel'] = channel

    year = int(args.year)
    default_config['year'] = year

    if year in (2017, 2018):
        args.weight_file = 'samples/xsec_2017.conf'

    if year == 2018:
        # FIXME: Need to update JEC when running on NanoAODv5
        default_config['jec'] = True

#     args.batch = True
    basename = os.path.basename(args.outputdir) + '_' + args.jet_type + '_' + channel + '_' + str(year)
    args.outputdir = os.path.join(os.path.dirname(args.outputdir), basename, 'data' if args.run_data else 'mc')
    args.jobdir = os.path.join('jobs_%s' % basename, 'data' if args.run_data else 'mc')

    if args.run_data:
        args.datasets = '%s/%s_%d_DATA.yaml' % (args.sample_dir, channel, year)
        args.extra_transfer = os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/JSON/%s' % golden_json[year])
        args.json = golden_json[year]
    else:
        args.datasets = '%s/%s_%d_MC.yaml' % (args.sample_dir, channel, year)

    if args.jet_type == 'ak15':
        args.cut = cut_dict_ak15[channel]
    else:
        args.cut = cut_dict_ak8[channel]

    args.imports = [('PhysicsTools.NanoHRTTools.producers.HeavyFlavSFTreeProducer', 'heavyFlavSFTreeFromConfig')]
    if not args.run_data:
        args.imports.extend([
            ('PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer', 'puAutoWeight_2017' if year == 2017 else 'puWeight_%d' % year),
            ])

    # data, or just nominal MC
    if args.run_data or not args.run_syst:
        cfg = copy.deepcopy(default_config)
        if args.run_data:
            cfg['data'] = True
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

#         # MET unclustEn up/down
#         for variation in ['up', 'down']:
#             syst_name = 'met_%s' % variation
#             logging.info('Start making %s trees...' % syst_name)
#             opts = copy.deepcopy(args)
#             cfg = copy.deepcopy(default_config)
#             cfg['met_unclustered'] = variation
#             opts.outputdir = os.path.join(os.path.dirname(opts.outputdir), syst_name)
#             opts.jobdir = os.path.join(os.path.dirname(opts.jobdir), syst_name)
#             run(opts, configs={hrt_cfgname: cfg})


def main():
    parser = get_arg_parser()

    parser.add_argument('--jet-type',
                        choices=['ak8', 'ak15'],
                        required=True,
                        help='Jet type: ak8, ak15'
                        )

    parser.add_argument('--channel',
                        choices=['photon', 'qcd', 'signal'],
                        required=True,
                        help='Channel: photon, qcd, signal'
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
                        help='Year: 2016, 2017, 2018, or comma separated list e.g., `2016,2017,2018`'
                        )

    parser.add_argument('--sample-dir',
                        type=str,
                        default='samples',
                        help='Directory of the sample list files. Default: %(default)s'
                        )

    args = parser.parse_args()

    if not (args.post or args.add_weight or args.merge):
        tar_cmssw()

    if ',' in args.year:
        years = [int(y) for y in args.year.split(',') if y]
        for year in years:
            for cat in ['data', 'mc']:
                opts = copy.deepcopy(args)
                if cat == 'data':
                    opts.run_data = True
                    opts.nfiles_per_job *= 2
                opts.inputdir = os.path.join(opts.inputdir.replace('_YEAR_', str(year)), cat)
                opts.year = year
                print(opts.inputdir, opts.year, opts.channel, 'data' if opts.run_data else 'mc', 'syst' if opts.run_syst else '')
                _process(opts)
    else:
        _process(args)


if __name__ == '__main__':
    main()

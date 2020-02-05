#!/usr/bin/env python
from __future__ import print_function

import os
import copy
from runPostProcessing import get_arg_parser, run, tar_cmssw
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')

hrt_cfgname = 'hrtSFTree_cfg.json'
default_config = {'data': False, 'channel': None, 'year': None, 'jec': False, 'jes': None, 'jes_source': '', 'jer': 'nominal', 'jmr': None, 'met_unclustered': None}
cut_dict = {
    'muon': 'Sum$(Muon_pt>55 && abs(Muon_eta)<2.4 && Muon_tightId && Muon_miniPFRelIso_all<0.10)>0 && nFatJet>0 && Sum$(abs(Jet_eta)<2.4 && Jet_btagDeepB>{DeepCSV_WP_M})>0',
    'photon': 'Sum$(Photon_pt>200)>0 && nFatJet>0',
    'qcd': 'Sum$((Jet_pt>25 && abs(Jet_eta)<2.4 && (Jet_jetId & 2)) * Jet_pt)>800 && nFatJet>0',
    }

golden_json = {
    2016: 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt',
    2017: 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt',
    2018: 'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt',
    }

def main():
    parser = get_arg_parser()

    parser.add_argument('--channel',
                        choices=['muon', 'photon', 'qcd'],
                        required=True,
                        help='Channel: muon, photon, qcd'
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
                        type=int,
                        choices=[2016, 2017, 2018],
                        required=True,
                        help='Year: 2016, 2017, 2018'
                        )

    args = parser.parse_args()

    channel = args.channel
    default_config['channel'] = channel

    year = args.year
    default_config['year'] = year

    if year in (2017, 2018):
        args.weight_file = 'samples/xsec_2017.conf'

    if year == 2018:
        # FIXME: Need to update JEC when running on NanoAODv5
        default_config['jec'] = True

    year_dep_cuts = {'DeepCSV_WP_M': {2016:0.6321, 2017:0.4941, 2018:0.4184}[year]}

    if not (args.post or args.add_weight or args.merge):
        tar_cmssw()

#     args.batch = True
    basename = os.path.basename(args.outputdir) + '_' + channel + '_' + str(year)
    args.outputdir = os.path.join(os.path.dirname(args.outputdir), basename, 'data' if args.run_data else 'mc')
    args.jobdir = os.path.join('jobs_%s' % basename, 'data' if args.run_data else 'mc')

    if args.run_data:
        args.datasets = 'samples/%s_%d_DATA.yaml' % (channel, year)
        args.json = os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/JSON/%s' % golden_json[year])
    else:
        args.datasets = 'samples/%s_%d_MC.yaml' % (channel, year)

    args.cut = cut_dict[channel].format(**year_dep_cuts)

    args.imports = [('PhysicsTools.NanoHRTTools.producers.hrtSFTreeProducer', 'hrtSFTreeFromConfig')]
    if not args.run_data:
        args.imports.extend([
            ('PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer', 'puAutoWeight_2017' if year == 2017 else 'puWeight_%d' % year),
            ])
        if args.channel == 'muon':
            args.imports.append(('PhysicsTools.NanoHRTTools.producers.topPtWeightProducer', 'topPtWeight')) #FIXME: year dependence?

    # data, or just nominal MC
    if args.run_data or not args.run_syst:
        cfg = copy.deepcopy(default_config)
        if args.run_data:
            cfg['data'] = True
        run(args, configs={hrt_cfgname: cfg})
        return

    # MC for syst.
    if args.run_syst:

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
        for variation in ['up', 'down']:
            syst_name = 'met_%s' % variation
            logging.info('Start making %s trees...' % syst_name)
            opts = copy.deepcopy(args)
            cfg = copy.deepcopy(default_config)
            cfg['met_unclustered'] = variation
            opts.outputdir = os.path.join(os.path.dirname(opts.outputdir), syst_name)
            opts.jobdir = os.path.join(os.path.dirname(opts.jobdir), syst_name)
            run(opts, configs={hrt_cfgname: cfg})


if __name__ == '__main__':
    main()

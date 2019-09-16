#!/usr/bin/env python
from __future__ import print_function

import os
import copy
from runPostProcessing import get_arg_parser, run, tar_cmssw
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')

hrt_cfgname = 'hrtSFTree_cfg.json'
default_config = {'data':False, 'channel': None, 'year': None, 'jec': False, 'jes': None, 'jes_source': '', 'jer': 'nominal', 'jmr':None, 'met_unclustered': None}

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

    if (year==2016):
	cut_dict = {
	    'muon': 'Sum$(Muon_pt>55 && Muon_tightId)>0 && (nCustomAK8Puppi+nCA15Puppi+nHOTVRPuppi)>0',	  
	    'photon': 'Sum$(Photon_pt>200)>0 && (nCustomAK8Puppi+nCA15Puppi+nHOTVRPuppi)>0',
	    'qcd': 'Sum$((Jet_pt>25 && abs(Jet_eta)<2.4 && (Jet_jetId & 2)) * Jet_pt)>800 && (nCustomAK8Puppi+nCA15Puppi+nHOTVRPuppi)>0',	   
    	}
	    
    elif(year==2017 or year==2018):
	cut_dict = {
    	'muon': 'Sum$(Muon_pt>55 && Muon_tightId)>0 && (nCustomAK8Puppi)>0', 
    	'photon': 'Sum$(Photon_pt>200)>0 && (nCustomAK8Puppi)>0',
    	'qcd': 'Sum$((Jet_pt>25 && abs(Jet_eta)<2.4 && (Jet_jetId & 2)) * Jet_pt)>800 && (nCustomAK8Puppi)>0',
   	 }
  

    if not (args.post or args.add_weight or args.merge):
        tar_cmssw()

#     args.batch = True
    basename = os.path.basename(args.outputdir) + '_' + channel + '_' + str(year)
    args.outputdir = os.path.join(os.path.dirname(args.outputdir), basename, 'data' if args.run_data else 'mc')
    args.jobdir = os.path.join('jobs_%s' % basename, 'data' if args.run_data else 'mc')
    args.datasets = 'samples/%s_%d.conf' % (channel, year)
    print (args.datasets)
    args.cut = cut_dict[channel]

    args.imports = [('PhysicsTools.NanoHRTTools.producers.hrtSFTreeProducer', 'hrtSFTreeFromConfig')]
    if not args.run_data:
        args.imports.extend([
            ('PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer', 'puAutoWeight_' + str(year)),
            ])
        if args.channel == 'muon':
            args.imports.append(('PhysicsTools.NanoHRTTools.producers.topPtWeightProducer', 'topPtWeight'))

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
#        opts.datasets = 'samples/%s_syst.conf' % channel
#        if not os.path.exists(opts.datasets):
#            opts.datasets = args.datasets
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

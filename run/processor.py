#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import time
import json
import argparse
import subprocess
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


def xrd(filepath):
    prefix = ''
    if filepath.startswith('/eos/cms'):
        prefix = 'root://eoscms.cern.ch/'
    elif filepath.startswith('/eos/uscms'):
        prefix = 'root://cmseos.fnal.gov/'
    if prefix:
        return prefix + '/' + filepath
    else:
        return filepath


def outputName(md, jobid):
    info = md['jobs'][jobid]
    return '{samp}_{idx}_tree.root'.format(samp=info['samp'], idx=info['idx'])


def main(args):

    # load job metadata
    with open(args.metadata) as fp:
        md = json.load(fp)

    # load modules
    modules = []
    for mod, names in md['imports']:
        import_module(mod)
        obj = sys.modules[mod]
        selnames = names.split(",")
        for name in dir(obj):
            if name[0] == "_":
                continue
            if name in selnames:
                print("Loading %s from %s " % (name, mod))
                modules.append(getattr(obj, name)())

    # run postprocessor
    outputname = outputName(md, args.jobid)
    p = PostProcessor(outputDir='.',
                      inputFiles=[xrd(f) for f in md['jobs'][args.jobid]['inputfiles']],
                      cut=md.get('cut'),
                      branchsel=os.path.basename(md['branchsel_in']),
                      modules=modules,
                      compression=md.get('compression', 'LZMA:9'),
                      friend=md.get('friend', False),
                      postfix=md.get('postfix'),
                      jsonInput=md.get('json'),
                      provenance=md.get('provenance', False),
                      haddFileName=None,
                      outputbranchsel=os.path.basename(md['branchsel_out']),
                      )
    p.run()

    # hadd files
    p = subprocess.Popen('haddnano.py %s *.root' % outputname, shell=True)
    p.communicate()
    if p.returncode != 0:
        raise RuntimeError('Hadd failed!')

    # keep only the hadd file
    for f in os.listdir('.'):
        if f.endswith('.root') and f != outputname:
            os.remove(f)

    # stage out
    if md['outputdir'].startswith('/eos'):
        cmd = 'xrdcp -np {outputname} {outputdir}/{outputname}'.format(outputname=outputname, outputdir=xrd(md['outputdir']))
        success = False
        for count in range(args.max_retry):
            p = subprocess.Popen(cmd, shell=True)
            p.communicate()
            if p.returncode == 0:
                success = True
                break
            else:
                time.sleep(args.sleep)
        if not success:
            raise RuntimeError("Stage out FAILED!")

        # clean up
        os.remove(outputname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NanoAOD postprocessing.')
    parser.add_argument('-m', '--metadata',
        default='metadata.json',
        help='Path to the metadata file. Default:%(default)s')
    parser.add_argument('--max-retry',
        type=int, default=3,
        help='Max number of retry for stageout. Default: %(default)s'
        )
    parser.add_argument('--sleep',
        type=int, default=120,
        help='Seconds to wait before retry stageout. Default: %(default)s'
        )
    parser.add_argument('jobid', type=int, help='Index of the output job.')

    args = parser.parse_args()
    main(args)

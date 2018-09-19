#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import json
import argparse
import re
import multiprocessing
import functools

import logging
logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s] %(levelname)s: %(message)s')


def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def tar_cmssw():
    cmsswdir = os.environ['CMSSW_BASE']
    cmsswtar = os.path.abspath(os.path.expandvars('$CMSSW_BASE/../CMSSW.tar.gz'))
    if os.path.exists(cmsswtar):
        ans = raw_input('CMSSW tarball %s already exists, remove? [yn] ' % cmsswtar)
        if ans.lower()[0] == 'y':
            os.remove(cmsswtar)
        else:
            return

    def exclude(tarinfo):
        exclude_patterns = ['/.git/', '/tmp/', '/jobs.*/', '/logs/', ]
        for pattern in exclude_patterns:
            if re.search(pattern, tarinfo.name):
                logging.debug('ignoring %s in the tarball', tarinfo.name)
                tarinfo = None
                break
        return tarinfo

    import tarfile
    with tarfile.open(cmsswtar, "w:gz") as tar:
        tar.add(cmsswdir, arcname=os.path.basename(cmsswdir), filter=exclude)


def create_metadata(args):
    '''
    create metadata

    Metadata is a dict including:
        - options
        - 'samples': (list)
        - 'inputfiles': (dict, sample -> files)
        - 'jobs': (list of dict)
            - jobitem: (dict, keys: 'samp', 'idx', 'inputfiles')
    '''

    arg_blacklist = ['metadata', 'select', 'ignore', 'site']
    md = {k: args.__dict__[k] for k in args.__dict__ if k not in arg_blacklist}

    md['samples'] = []
    md['inputfiles'] = {}
    md['jobs'] = []

    # discover all the datasets
    for samp in os.listdir(args.inputdir):
        if args.select:
            sels = args.select.split(',')
            match = False
            for s in sels:
                if re.search(s, samp):
                    logging.debug('Selecting dataset %s', samp)
                    match = True
                    break
            if not match:
                continue
        elif args.ignore:
            vetoes = args.ignore.split(',')
            match = False
            for v in vetoes:
                if re.search(v, samp):
                    logging.debug('Ignoring dataset %s', samp)
                    match = True
                    break
            if match:
                continue
        md['samples'].append(samp)

    # sort the samples
    md['samples'] = natural_sort(md['samples'])

    # discover the files
    for samp in md['samples']:
        # get all files of this sample
        md['inputfiles'][samp] = []
        sampdir = os.path.join(args.inputdir, samp)
        for dp, dn, filenames in os.walk(sampdir):
            if 'failed' in dp:
                continue
            for f in filenames:
                if not f.endswith('.root'):
                    continue
                md['inputfiles'][samp].append(os.path.join(dp, f))

        # sort the input list
        md['inputfiles'][samp] = natural_sort(md['inputfiles'][samp])

        # create jobs
        for idx, chunk in enumerate(get_chunks(md['inputfiles'][samp], args.nfiles_per_job)):
            md['jobs'].append({'samp': samp, 'idx': idx, 'inputfiles': chunk})

    return md


def load_metadata(args):
    metadatafile = os.path.join(args.jobdir, args.metadata)
    with open(metadatafile) as f:
        md = json.load(f)
    return md


def submit(args):
    scriptfile = os.path.join(os.path.dirname(__file__), 'run_postproc_condor.sh')
    macrofile = os.path.join(os.path.dirname(__file__), 'processor.py')
    metadatafile = os.path.join(args.jobdir, args.metadata)
    joboutputdir = os.path.join(args.outputdir, 'pieces')

    if not args.resubmit:
        # create jobdir
        if os.path.exists(args.jobdir):
            ans = raw_input('jobdir %s already exists, remove? [yn] ' % args.jobdir)
            if ans.lower()[0] == 'y':
                import shutil
                shutil.rmtree(args.jobdir)
            else:
                sys.exit(1)
        os.makedirs(args.jobdir)

        # create outputdir
        if not os.path.exists(joboutputdir):
            os.makedirs(joboutputdir)

        # create metadata file
        md = create_metadata(args)
        with open(metadatafile, 'w') as f:
            json.dump(md, f, ensure_ascii=True, indent=2, sort_keys=True)

        # create CMSSW tarball
        tar_cmssw()

        njobs = len(md['jobs'])
        jobids = [str(jobid) for jobid in range(njobs)]
        jobids_file = os.path.join(args.jobdir, 'submit.txt')

    else:
        # resubmit
        jobids = []
        jobids_file = os.path.join(args.jobdir, 'resubmit.txt')
        log_files = [f for f in os.listdir(args.jobdir) if f.endswith('.log')]
        for fn in log_files:
            with open(os.path.join(args.jobdir, fn)) as logfile:
                errormsg = None
                for line in reversed(logfile.readlines()):
                    if 'Job removed' in line or 'aborted' in line:
                        errormsg = line
                    if 'Job submitted from host' in line:
                        # if seeing this first: the job has been resubmited
                        break
                    if 'return value' in line:
                        if 'return value 0' not in line:
                            errormsg = line
                        break
                if errormsg:
                    logging.debug(fn + '\n   ' + errormsg)
                    jobids.append(fn.split('.')[0])
                    assert jobids[-1].isdigit()

    with open(jobids_file, 'w') as f:
        f.write('\n'.join(jobids))

    condordesc = '''\
universe              = vanilla
requirements          = (Arch == "X86_64") && (OpSys == "LINUX")
request_disk          = 10000000
executable            = {scriptfile}
arguments             = $(jobid)
transfer_input_files  = {macro},{cmsswtar},{metadatafile}
output                = {jobdir}/$(jobid).out
error                 = {jobdir}/$(jobid).err
log                   = {jobdir}/$(jobid).log
use_x509userproxy     = true
Should_Transfer_Files = YES
initialdir            = {outputdir}
WhenToTransferOutput  = ON_EXIT
want_graceful_removal = true
on_exit_remove        = (ExitBySignal == False) && (ExitCode == 0)
on_exit_hold          = ( (ExitBySignal == True) || (ExitCode != 0) )
on_exit_hold_reason   = strcat("Job held by ON_EXIT_HOLD due to ", ifThenElse((ExitBySignal == True), "exit by signal", strcat("exit code ",ExitCode)), ".")
periodic_release      = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 10*60)
{site}

queue jobid from {jobids_file}
'''.format(scriptfile=os.path.abspath(scriptfile),
           macro=os.path.abspath(macrofile),
           metadatafile=os.path.abspath(metadatafile),
           cmsswtar=os.path.abspath(os.path.expandvars('$CMSSW_BASE/../CMSSW.tar.gz')),
           jobdir=os.path.abspath(args.jobdir),
           outputdir=joboutputdir,
           jobids_file=os.path.abspath(jobids_file),
           site='+DESIRED_Sites = "%s"' % args.site if args.site else '',
    )
    condorfile = os.path.join(args.jobdir, 'submit.cmd')
    with open(condorfile, 'w') as f:
        f.write(condordesc)

    print('Run the following command to submit the jobs:\n  condor_submit {condorfile}'.format(condorfile=condorfile))


def run_merge(args):
    import subprocess
    md = load_metadata(args)
    for samp in md['samples']:
        cmd = 'haddnano.py {outdir}/{samp}_tree.root {outdir}/pieces/{samp}_*_tree.root'.format(outdir=args.outputdir, samp=samp)
        logging.debug('...' + cmd)
        p = subprocess.Popen(cmd, shell=True)
        p.communicate()
        if p.returncode != 0:
            raise RuntimeError('Hadd failed on %s!' % samp)


def run_all(args):
    pass
#     md, njobs = update_metadata(args)
#     pool = multiprocessing.Pool(args.nproc)
#     pool.map(
#         functools.partial(writeData, md, args.outputdir, batch_mode=False,
#                           test_sample=args.test_sample, events=args.events_per_file, dryrun=args.dryrun),
#         range(njobs)
#         )


def main():
    parser = argparse.ArgumentParser('Preprocess ntuples')
    parser.add_argument('inputdir',
        help='Input diretory.'
    )
    parser.add_argument('outputdir',
        help='Output directory'
    )
    parser.add_argument('-m', '--metadata',
        default='metadata.json',
        help='Metadata json file. Default: %(default)s'
    )
    parser.add_argument('-t', '--submittype',
        default='condor', choices=['interactive', 'condor'],
        help='Method of job submission. [Default: %(default)s]'
        )
    parser.add_argument('--resubmit',
        action='store_true', default=False,
        help='Resubmit failed jobs. Default: %(default)s'
    )
    parser.add_argument('-j', '--jobdir',
        default='jobs',
        help='Directory for job files. [Default: %(default)s]'
        )
    parser.add_argument('--select',
        default='',
        help='Selected datasets, common separated regex. [Default: %(default)s]'
    )
    parser.add_argument('--ignore',
        default='',
        help='Ignored datasets, common separated regex. [Default: %(default)s]'
    )
#     parser.add_argument('--nproc',
#         type=int, default=8,
#         help='Number of jobs to run in parallel. Default: %(default)s'
#     )
    parser.add_argument('-n', '--nfiles-per-job',
        type=int, default=10,
        help='Number of input files to process in one job. Default: %(default)s'
    )
    parser.add_argument('--dryrun',
        action='store_true', default=False,
        help='Do not convert -- only produce metadata. Default: %(default)s'
    )
    parser.add_argument('--site',
        default='',
        help='Specify sites for condor submission. Default: %(default)s'
    )

    parser.add_argument("-s", "--postfix", dest="postfix", default=None, help="Postfix which will be appended to the file name (default: _Friend for friends, _Skim for skims)")
    parser.add_argument("-J", "--json", dest="json", default=None, help="Select events using this JSON file")
    parser.add_argument("-c", "--cut", dest="cut", default=None, help="Cut string")
    parser.add_argument("-b", "--branch-selection", dest="branchsel", default=None, help="Branch selection")
    parser.add_argument("--friend", dest="friend", action="store_true", default=False, help="Produce friend trees in output (current default is to produce full trees)")
    parser.add_argument("-I", "--import", dest="imports", default=[], action="append", nargs=2, help="Import modules (python package, comma-separated list of ")
    parser.add_argument("-z", "--compression", dest="compression", default=("LZMA:9"), help="Compression: none, or (algo):(level) ")

    parser.add_argument('--merge',
        action='store_true', default=False,
        help='Merge output files. Default: %(default)s'
    )

    args = parser.parse_args()
#     print(args)

    if args.merge:
        run_merge(args)

    if args.submittype == 'interactive':
        run_all(args)
    elif args.submittype == 'condor':
        submit(args)


if __name__ == '__main__':
    main()

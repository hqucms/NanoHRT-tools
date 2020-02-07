#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import json
import re
import shutil

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


def get_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def sname(sample_or_dataset_name):
    if sample_or_dataset_name[0] == '/':
        return sample_or_dataset_name.split('/')[1]
    else:
        return sample_or_dataset_name


def get_filenames(dataset, retry=3):
    """Return files for given DAS query via dasgoclient"""
    import subprocess
    import time
    query = 'file dataset=%s' % dataset
    if dataset.endswith('/USER'):
        query += ' instance=prod/phys03'
    cmd = ['dasgoclient', '-query', query, '-json']
    retry_count = 0
    while True:
        logging.info('Querying DAS:\n  %s' % ' '.join(cmd) + '' if retry_count == 0 else '\n... retry %d ...' % retry_count)
        if retry_count > 0:
            time.sleep(3)
        retry_count += 1
        if retry_count > retry:
            raise RuntimeError('Failed to retrieve file names from DAS for: %s' % dataset)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outs, errs = proc.communicate()
        if errs:
            logging.error('DAS error: %s' % errs)
            continue
        else:
            files = []
            for row in json.loads(outs):
                for rec in row.get('file', []):
                    fname = rec.get('name', '')
                    if fname:
                        files.append(str(fname))
            logging.info('Found %d files for %s' % (len(files), dataset))
            return files


def add_weight_branch(file, xsec, lumi=1000., treename='Events', wgtbranch='xsecWeight'):
    from array import array
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True

    def _get_sum(tree, wgtvar):
        htmp = ROOT.TH1D('htmp', 'htmp', 1, 0, 10)
        tree.Project('htmp', '1.0', wgtvar)
        return float(htmp.Integral())

    def _fill_const_branch(tree, branch_name, buff, lenVar=None):
        if lenVar is not None:
            b = tree.Branch(branch_name, buff, '%s[%s]/F' % (branch_name, lenVar))
            b_lenVar = tree.GetBranch(lenVar)
            buff_lenVar = array('I', [0])
            b_lenVar.SetAddress(buff_lenVar)
        else:
            b = tree.Branch(branch_name, buff, branch_name + '/F')

        b.SetBasketSize(tree.GetEntries() * 2)  # be sure we do not trigger flushing
        for i in xrange(tree.GetEntries()):
            if lenVar is not None:
                b_lenVar.GetEntry(i)
            b.Fill()

        b.ResetAddress()
        if lenVar is not None:
            b_lenVar.ResetAddress()

    f = ROOT.TFile(file, 'UPDATE')
    run_tree = f.Get('Runs')
    tree = f.Get(treename)

    # fill cross section weights to the 'Events' tree
    sumwgts = _get_sum(run_tree, 'genEventSumw')
    xsecwgt = xsec * lumi / sumwgts
    xsec_buff = array('f', [xsecwgt])
    _fill_const_branch(tree, wgtbranch, xsec_buff)

    # fill LHE weight re-normalization factors
    if tree.GetBranch('LHEScaleWeight'):
        run_tree.GetEntry(0)
        nScaleWeights = run_tree.nLHEScaleSumw
        scale_weight_norm_buff = array('f', [sumwgts / _get_sum(run_tree, 'LHEScaleSumw[%d]*genEventSumw' % i) for i in xrange(nScaleWeights)])
        logging.info('LHEScaleWeightNorm: ' + str(scale_weight_norm_buff))
        _fill_const_branch(tree, 'LHEScaleWeightNorm', scale_weight_norm_buff, lenVar='nLHEScaleWeight')

    if tree.GetBranch('LHEPdfWeight'):
        run_tree.GetEntry(0)
        nPdfWeights = run_tree.nLHEPdfSumw
        pdf_weight_norm_buff = array('f', [sumwgts / _get_sum(run_tree, 'LHEPdfSumw[%d]*genEventSumw' % i) for i in xrange(nPdfWeights)])
        logging.info('LHEPdfWeightNorm: ' + str(pdf_weight_norm_buff))
        _fill_const_branch(tree, 'LHEPdfWeightNorm', pdf_weight_norm_buff, lenVar='nLHEPdfWeight')

    tree.Write(treename, ROOT.TObject.kOverwrite)
    f.Close()


def load_dataset_file(dataset_file):
    import yaml
    with open(dataset_file) as f:
        d = yaml.safe_load(f)

    outtree_to_samples = {}
    samp_to_datasets = {}
    for outtree_name in d:
        outtree_to_samples[outtree_name] = []
        for samp_or_list in d[outtree_name]:
            if isinstance(samp_or_list, list):
                samp = sname(samp_or_list[0])
                sample_list = samp_or_list
            else:
                samp = sname(samp_or_list)
                sample_list = [samp_or_list]
            outtree_to_samples[outtree_name].append(samp)
            samp_to_datasets[samp] = sample_list
    return outtree_to_samples, samp_to_datasets


def parse_sample_xsec(cfgfile):
    xsec_dict = {}
    with open(cfgfile) as f:
        for l in f:
            l = l.strip()
            if not l or l.startswith('#'):
                continue
            pieces = l.split()
            samp = None
            xsec = None
            isData = False
            for s in pieces:
                if '/MINIAOD' in s or '/NANOAOD' in s:
                    samp = s.split('/')[1]
                    if '/MINIAODSIM' not in s and '/NANOAODSIM' not in s:
                        isData = True
                        break
                else:
                    try:
                        xsec = float(s)
                    except ValueError:
                        try:
                            import numexpr
                            xsec = numexpr.evaluate(s).item()
                        except:
                            pass
            if samp is None:
                logging.warning('Ignore line:\n%s' % l)
            elif not isData and xsec is None:
                logging.error('Cannot find cross section:\n%s' % l)
            else:
                if samp in xsec_dict and xsec_dict[samp] != xsec:
                    raise RuntimeError('Inconsistent entries for sample %s' % samp)
                xsec_dict[samp] = xsec
                if 'PSweights_' in samp:
                    xsec_dict[samp.replace('PSweights_', '')] = xsec
    return xsec_dict


def tar_cmssw(batchMode=False):
    cmsswdir = os.environ['CMSSW_BASE']
    cmsswtar = os.path.abspath(os.path.expandvars('$CMSSW_BASE/../CMSSW.tar.gz'))
    if os.path.exists(cmsswtar):
        if batchMode:
            return
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

    arg_blacklist = ['metadata', 'select', 'ignore', 'site', 'datasets']
    md = {k: args.__dict__[k] for k in args.__dict__ if k not in arg_blacklist}

    md['samples'] = []
    md['inputfiles'] = {}
    md['jobs'] = []

    def select_sample(dataset):
        samp = sname(dataset)
        keep = True
        if args.select:
            sels = args.select.split(',')
            match = False
            for s in sels:
                if re.search(s, samp):
                    logging.debug('Selecting dataset %s', dataset)
                    match = True
                    break
            if not match:
                keep = False
        elif args.ignore:
            vetoes = args.ignore.split(',')
            match = False
            for v in vetoes:
                if re.search(v, samp):
                    logging.debug('Ignoring dataset %s', dataset)
                    match = True
                    break
            if match:
                keep = False
        return keep

    _, samp_to_datasets = load_dataset_file(args.datasets)

    # discover all the datasets
    if args.inputdir:
        found_samples = os.listdir(args.inputdir)
        for samp in samp_to_datasets:
            filelist = []
            for dataset in samp_to_datasets[samp]:
                if sname(dataset) not in found_samples:
                    logging.warning('Cannot find dataset %s in the input dir %s' % (dataset, args.inputdir))
                    continue
                if not select_sample(dataset):
                    continue
                sampdir = os.path.join(args.inputdir, sname(dataset))
                for dp, dn, filenames in os.walk(sampdir):
                    if 'failed' in dp:
                        continue
                    for f in filenames:
                        if not f.endswith('.root'):
                            continue
                        if os.path.getsize(os.path.join(dp, f)) < 1000:
                            if '-' not in samp and '_' not in samp:
                                raise RuntimeError('Empty data file %s' % os.path.join(dp, f))
                            else:
                                logging.warning('Ignoring empty MC file %s' % os.path.join(dp, f))
                                continue
                        filelist.append(os.path.join(dp, f))
            if len(filelist):
                md['samples'].append(samp)
                md['inputfiles'][samp] = filelist
    else:
        # use remote files
        for samp in samp_to_datasets:
            filelist = []
            dataset0 = None
            for dataset in samp_to_datasets[samp]:
                if dataset0 is None:
                    dataset0 = dataset.split('/')[1]
                else:
                    if dataset0 != dataset.split('/')[1]:
                        raise RuntimeError('Inconsistent dataset for samp `%s`: `%s` vs `%s`' % (samp, dataset0, dataset))
                if select_sample(dataset):
                    filelist.extend(get_filenames(dataset))
            if len(filelist):
                md['samples'].append(samp)
                md['inputfiles'][samp] = filelist

    # sort the samples
    md['samples'] = natural_sort(md['samples'])

    # discover the files
    for samp in md['samples']:
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


def check_job_status(args):
    metadatafile = os.path.join(args.jobdir, args.metadata)
    with open(metadatafile) as f:
        md = json.load(f)
    njobs = len(md['jobs'])
    jobids = {'running': [], 'failed': [], 'completed': []}
    for jobid in range(njobs):
        logpath = os.path.join(args.jobdir, '%d.log' % jobid)
        if not os.path.exists(logpath):
            logging.debug('Cannot find log file %s' % logpath)
            jobids['failed'].append(str(jobid))
            continue
        with open(logpath) as logfile:
            errormsg = None
            finished = False
            for line in reversed(logfile.readlines()):
                if 'Job removed' in line or 'aborted' in line:
                    errormsg = line
                if 'Job submitted from host' in line:
                    # if seeing this first: the job has been resubmited
                    break
                if 'return value' in line:
                    if 'return value 0' in line:
                        finished = True
                    else:
                        errormsg = line
                    break
            if errormsg:
                logging.debug(logpath + '\n   ' + errormsg)
                jobids['failed'].append(str(jobid))
            else:
                if finished:
                    jobids['completed'].append(str(jobid))
                else:
                    jobids['running'].append(str(jobid))
    assert sum(len(jobids[k]) for k in jobids) == njobs
    all_completed = len(jobids['completed']) == njobs
    info = {k: len(jobids[k]) for k in jobids if len(jobids[k])}
    logging.info('Job %s status: ' % args.jobdir + str(info))
    return all_completed, jobids


def submit(args, configs):
    logging.info('Preparing jobs...\n  - modules: %s\n  - cut: %s\n  - outputdir: %s' % (str(args.imports), args.cut, args.outputdir))

    scriptfile = os.path.join(os.path.dirname(__file__), 'run_postproc_condor.sh')
    macrofile = os.path.join(os.path.dirname(__file__), 'processor.py')
    metadatafile = os.path.join(args.jobdir, args.metadata)
    joboutputdir = os.path.join(args.outputdir, 'pieces')

    # create config file for the scripts
    configfiles = []
    if configs is not None:
        for cfgname in configs:
            cfgpath = os.path.join(args.jobdir, cfgname)
            configfiles.append(cfgpath)

    if not args.resubmit:
        # create jobdir
        if os.path.exists(args.jobdir):
            if args.batch:
                logging.warning('jobdir %s already exists! Will not submit new jobs!' % args.jobdir)
                return
            ans = raw_input('jobdir %s already exists, remove? [yn] ' % args.jobdir)
            if ans.lower()[0] == 'y':
                shutil.rmtree(args.jobdir)
            else:
                sys.exit(1)
        os.makedirs(args.jobdir)

        # create outputdir
        if os.path.exists(joboutputdir):
            if not args.batch:
                ans = raw_input('outputdir %s already exists, continue? [yn] ' % joboutputdir)
                if ans.lower()[0] == 'n':
                    sys.exit(1)
        else:
            os.makedirs(joboutputdir)

        # create config file for the scripts
        if configs is not None:
            for cfgname, cfgpath in zip(configs, configfiles):
                with open(cfgpath, 'w') as f:
                    json.dump(configs[cfgname], f, ensure_ascii=True, indent=2, sort_keys=True)
                shutil.copy2(cfgpath, joboutputdir)

        # create metadata file
        md = create_metadata(args)
        with open(metadatafile, 'w') as f:
            json.dump(md, f, ensure_ascii=True, indent=2, sort_keys=True)
        shutil.copy2(metadatafile, joboutputdir)

        # create CMSSW tarball
        tar_cmssw(args.batch)

        njobs = len(md['jobs'])
        jobids = [str(jobid) for jobid in range(njobs)]
        jobids_file = os.path.join(args.jobdir, 'submit.txt')

    else:
        # resubmit
        jobids = check_job_status(args)[1]['failed']
        jobids_file = os.path.join(args.jobdir, 'resubmit.txt')

    with open(jobids_file, 'w') as f:
        f.write('\n'.join(jobids))

    # prepare the list of files to transfer
    files_to_transfer = [os.path.expandvars('$CMSSW_BASE/../CMSSW.tar.gz'), macrofile, metadatafile] + configfiles
    if args.branchsel_in:
        files_to_transfer.append(args.branchsel_in)
    if args.branchsel_out:
        files_to_transfer.append(args.branchsel_out)
    if args.extra_transfer:
        files_to_transfer += args.extra_transfer.split(',')
    files_to_transfer = [os.path.abspath(f) for f in files_to_transfer]

    condordesc = '''\
universe              = vanilla
requirements          = (Arch == "X86_64") && (OpSys == "LINUX")
request_disk          = 10000000
executable            = {scriptfile}
arguments             = $(jobid)
transfer_input_files  = {files_to_transfer}
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
           files_to_transfer=','.join(files_to_transfer),
           jobdir=os.path.abspath(args.jobdir),
           outputdir=joboutputdir,
           jobids_file=os.path.abspath(jobids_file),
           site='+DESIRED_Sites = "%s"' % args.site if args.site else '',
    )
    condorfile = os.path.join(args.jobdir, 'submit.cmd')
    with open(condorfile, 'w') as f:
        f.write(condordesc)

    cmd = 'condor_submit {condorfile}'.format(condorfile=condorfile)
    print('Run the following command to submit the jobs:\n  %s' % cmd)
    if args.batch:
        import subprocess
        subprocess.Popen(cmd, shell=True).communicate()


def run_add_weight(args):
    if args.weight_file:
        xsec_dict = parse_sample_xsec(args.weight_file)

    import subprocess
    md = load_metadata(args)
    parts_dir = os.path.join(args.outputdir, 'parts')
    if not os.path.exists(parts_dir):
        os.makedirs(parts_dir)
    for samp in md['samples']:
        outfile = '{parts_dir}/{samp}_tree.root'.format(parts_dir=parts_dir, samp=samp)
        cmd = 'haddnano.py {outfile} {outputdir}/pieces/{samp}_*_tree.root'.format(outfile=outfile, outputdir=args.outputdir, samp=samp)
        logging.debug('...' + cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = p.communicate()[0]
        log_lower = log.lower()
        if 'error' in log_lower or 'fail' in log_lower:
            logging.error(log)
        if p.returncode != 0:
            raise RuntimeError('Hadd failed on %s!' % samp)

        # add weight
        if args.weight_file:
            try:
                xsec = xsec_dict[samp]
                if xsec is not None:
                    logging.info('Adding xsec weight to file %s, xsec=%f' % (outfile, xsec))
                    add_weight_branch(outfile, xsec)
            except KeyError as e:
                if '-' not in samp and '_' not in samp:
                    # data
                    logging.info('Not adding weight to sample %s' % samp)
                else:
                    raise e


def run_merge(args):
    import subprocess

    parts_dir = os.path.join(args.outputdir, 'parts')
    allfiles = [f for f in os.listdir(parts_dir) if f.endswith('.root')]
    merge_dict = {}  # outname: expected files
    merge_dict_found = {}  # outname: [infile list]
    outtree_to_samples, _ = load_dataset_file(args.datasets)

    for outtree_name in outtree_to_samples:
        outname = '%s_tree.root' % outtree_name
        merge_dict[outname] = []
        merge_dict_found[outname] = []
        for samp in outtree_to_samples[outtree_name]:
            fname = samp + '_tree.root'
            merge_dict[outname].append(os.path.join(parts_dir, fname))
            if fname in allfiles:
                merge_dict_found[outname].append(os.path.join(parts_dir, fname))

    for outname in merge_dict:
        if len(merge_dict_found[outname]) == 0:
            logging.warning('Ignore %s as no input files are found.' % outname)
            continue
        if len(merge_dict_found[outname]) != len(merge_dict[outname]):
            raise RuntimeError('Incomplete files for merging, missing: %s' % str(set(merge_dict[outname]) - set(merge_dict_found[outname])))

        if len(merge_dict_found[outname]) == 1:
            os.rename(list(merge_dict_found[outname])[0], os.path.join(args.outputdir, outname))
        else:
            cmd = 'haddnano.py {outfile} {infiles}'.format(outfile=os.path.join(args.outputdir, outname), infiles=' '.join(merge_dict_found[outname]))
            logging.debug('...' + cmd)
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log = p.communicate()[0]
            if p.returncode != 0:
                raise RuntimeError('Hadd failed on %s!' % outname)
            log_lower = log.lower()
            if 'error' in log_lower or 'fail' in log_lower:
                logging.error(log)


def run_all(args):
    import multiprocessing
    import functools
    pass
#     md, njobs = update_metadata(args)
#     pool = multiprocessing.Pool(args.nproc)
#     pool.map(
#         functools.partial(writeData, md, args.outputdir, batch_mode=False,
#                           test_sample=args.test_sample, events=args.events_per_file, dryrun=args.dryrun),
#         range(njobs)
#         )


def get_arg_parser():
    import argparse
    parser = argparse.ArgumentParser('Preprocess ntuples')
    parser.add_argument('-i', '--inputdir', default=None,
        help='Input diretory.'
    )
    parser.add_argument('-o', '--outputdir', required=True,
        help='Output directory'
    )
    parser.add_argument('-m', '--metadata',
        default='metadata.json',
        help='Metadata json file. Default: %(default)s'
    )
    parser.add_argument('--extra-transfer',
        default=None,
        help='Extra files to transfer, common separated list. Default: %(default)s'
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
    parser.add_argument('-d', '--datasets', required=False,
        default='',
        help='Path to the dataset list file. [Default: %(default)s]'
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
    parser.add_argument('--add-weight',
        action='store_true', default=False,
        help='Merge output files of the same dataset and add cross section weight using the file specified in --weight-file. Default: %(default)s'
    )
    parser.add_argument('-w', '--weight-file',
        default='samples/xsec.conf',
        help='File with xsec of each sample. If empty, xsec wgt will not be added. Default: %(default)s'
    )
    parser.add_argument('--merge',
        action='store_true', default=False,
        help='Merge output files of different sample as specified in --datasets file. Default: %(default)s'
    )
    parser.add_argument('--post',
        action='store_true', default=False,
        help='Add weight and merge. Default: %(default)s'
    )
    parser.add_argument('--batch',
        action='store_true', default=False,
        help='Batch mode, do not ask for confirmation and submit the jobs directly. Default: %(default)s'
    )

    # preserve the options in nano_postproc.py
    parser.add_argument("-s", "--postfix", dest="postfix", default=None, help="Postfix which will be appended to the file name (default: _Friend for friends, _Skim for skims)")
    parser.add_argument("-J", "--json", dest="json", default=None, help="Select events using this JSON file")
    parser.add_argument("-c", "--cut", dest="cut", default=None, help="Cut string")
    parser.add_argument("--bi", "--branch-selection-input", dest="branchsel_in", default='keep_and_drop_input.txt', help="Branch selection input")
    parser.add_argument("--bo", "--branch-selection-output", dest="branchsel_out", default='keep_and_drop_output.txt', help="Branch selection output")
    parser.add_argument("--friend", dest="friend", action="store_true", default=False, help="Produce friend trees in output (current default is to produce full trees)")
    parser.add_argument("-I", "--import", dest="imports", default=[], action="append", nargs=2, help="Import modules (python package, comma-separated list of ")
    parser.add_argument("-z", "--compression", dest="compression", default=("LZ4:4"), help="Compression: none, or (algo):(level) ")
    parser.add_argument("-P", "--prefetch", dest="prefetch", action="store_true", default=False, help="Prefetch input files locally instead of accessing them via xrootd")
    parser.add_argument("--long-term-cache", dest="longTermCache", action="store_true", default=False, help="Keep prefetched files across runs instead of deleting them at the end")
    parser.add_argument("-N", "--max-entries", dest="maxEntries", type=int, default=None, help="Maximum number of entries to process from any single given input tree")
    parser.add_argument("--first-entry", dest="firstEntry", type=int, default=0, help="First entry to process in the three (to be used together with --max-entries)")
    parser.add_argument("--justcount", dest="justcount", default=False, action="store_true", help="Just report the number of selected events")

    return parser


def run(args, configs=None):

    if args.post:
        args.add_weight = True
        args.merge = True

    if args.add_weight:
        all_completed, _ = check_job_status(args)
        if not all_completed:
            ans = raw_input('Warning! There are jobs failed or still running. Continue adding weights? [yn] ')
            if ans.lower()[0] != 'y':
                sys.exit()
        run_add_weight(args)

    if args.merge:
        run_merge(args)

    if args.add_weight or args.merge:
        return

    if args.submittype == 'interactive':
        run_all(args, configs)
    elif args.submittype == 'condor':
        submit(args, configs)


if __name__ == '__main__':
    parser = get_arg_parser()
    args = parser.parse_args()
#     print(args)

    run(args)

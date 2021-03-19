import os
import numpy as np
import awkward
import onnxruntime
import json
import pandas as pd
import traceback
import logging

from .makeInputs import ParticleNetTagInfoMaker


def configLogger(name, loglevel=logging.INFO, filename=None):
    # define a Handler which writes INFO messages or higher to the sys.stderr
    logger = logging.getLogger(name)
    logger.setLevel(loglevel)
    console = logging.StreamHandler()
    console.setLevel(loglevel)
    console.setFormatter(logging.Formatter(
        '[%(asctime)s] %(levelname)s: %(message)s'))
    logger.addHandler(console)
    if filename:
        logfile = logging.FileHandler(filename)
        logfile.setLevel(loglevel)
        logfile.setFormatter(logging.Formatter(
            '[%(asctime)s] %(levelname)s: %(message)s'))
        logger.addHandler(logfile)


def md5(fname):
    '''https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file'''
    import hashlib
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def copyFileEOS(source, destination, max_retry=1, sleep=10):
    import subprocess
    import time

    cmd = 'xrdcp --silent -p -f {source} {destination}'.format(source=source, destination=destination)
    print(cmd)
    success = False
    for count in range(max_retry):
        p = subprocess.Popen(cmd, shell=True)
        p.communicate()
        if p.returncode == 0:
            success = True
            break
        else:
            time.sleep(sleep)
    if not success:
        raise RuntimeError("FAILED to copy file %s!" % source)


def _pad(a, min_length, max_length, value=0, dtype='float32'):
    if len(a) > max_length:
        return a[:max_length].astype(dtype)
    elif len(a) < min_length:
        x = (np.ones(min_length) * value).astype(dtype)
        x[:len(a)] = a.astype(dtype)
        return x
    else:
        return a.astype('float32')


logger = logging.getLogger('NanoNN')
configLogger('NanoNN', loglevel=logging.INFO)


class ParticleNetJetTagsProducer(object):

    def __init__(self, model_path, preprocess_path, version, cache_suffix, debug=False):
        self.debug = debug
        model_path = model_path.format(version=version)
        preprocess_path = preprocess_path.format(version=version)
        with open(preprocess_path) as fp:
            self.prep_params = json.load(fp)
        logger.info('Loading model %s' % model_path)
        self.sess = onnxruntime.InferenceSession(model_path)
        self.ver = version
        self.cache_suffix = cache_suffix
        self.md5 = md5(model_path)

    def _preprocess(self, taginfo, eval_flags=None):
        data = {}
        counts = None
        for group_name in self.prep_params['input_names']:
            data[group_name] = []
            info = self.prep_params[group_name]
            for var in info['var_names']:
                a = taginfo[var].copy()
                if eval_flags is not None:
                    a = a[eval_flags]
                if counts is None:
                    counts = a.counts
                else:
                    assert(np.array_equal(counts, a.counts))
                a = (a - info['var_infos'][var]['median']) * info['var_infos'][var]['norm_factor']
                a = a.flatten().pad(info['var_length'], clip=True).fillna(0).regular()
                a = np.clip(a, info['var_infos'][var].get('lower_bound', -5),
                            info['var_infos'][var].get('upper_bound', 5))
                if self.debug:
                    print(var, a)
                data[group_name].append(a.astype('float32'))
            data[group_name] = np.nan_to_num(np.stack(data[group_name], axis=1))
        return data, counts

    def predict(self, taginfo, eval_flags=None):
        data, counts = self._preprocess(taginfo, eval_flags)
        preds = self.sess.run([], data)[0]
        outputs = {
            flav: awkward.JaggedArray.fromcounts(counts, preds[:, i]) for i,
            flav in enumerate(self.prep_params['output_names'])}
        return outputs

    def predict_one(self, taginfo, entry_idx, jet_idx, jet=None):
        data = {}
        for group_name in self.prep_params['input_names']:
            data[group_name] = []
            info = self.prep_params[group_name]
            for var in info['var_names']:
                a = taginfo[var][entry_idx][jet_idx].copy()
                a = (a - info['var_infos'][var]['median']) * info['var_infos'][var]['norm_factor']
                a = np.clip(a, info['var_infos'][var].get('lower_bound', -5),
                            info['var_infos'][var].get('upper_bound', 5))
                try:
                    a = _pad(a, min_length=info['min_length'], max_length=info['max_length'])
                except KeyError:
                    a = _pad(a, min_length=info['var_length'], max_length=info['var_length'])
                if self.debug:
                    print(var, a)
                data[group_name].append(a.astype('float32'))
            data[group_name] = np.nan_to_num(np.expand_dims(np.stack(data[group_name], axis=0), 0))
        preds = self.sess.run([], data)[0]
        outputs = {flav: preds[0, i] for i, flav in enumerate(self.prep_params['output_names'])}
        if self.debug:
            p4 = taginfo['_jetp4'][entry_idx][jet_idx]
            print('pt,eta,phi', (jet.pt, jet.eta, jet.phi), (p4.pt, p4.eta, p4.phi))
            print('outputs', outputs)
        return outputs

    def load_cache(self, inputFile):
        self.cache_fullpath = inputFile.GetName().replace('.root', '.%s%s.h5' % (self.cache_suffix, self.ver))
        self.cachefile = os.path.basename(self.cache_fullpath)
        try:
            copyFileEOS(self.cache_fullpath, self.cachefile)
            self._cache_df = pd.read_hdf(self.cachefile, key=self.md5)
            self._cache_dict = self._cache_df.set_index(['event', 'jetidx']).to_dict(orient='index')
            logger.info('Loaded cache from %s' % self.cache_fullpath)
        except KeyError:
            raise
        except Exception:
            logger.warning('Cannot load the cache -- Will run the model from scratch...')
            self._cache_df = None
        self._cache_list = []
        return self._cache_df is not None

    def update_cache(self):
        if len(self._cache_list) > 0:
            logger.info('Updating the cache file to include %d new entries...' % len(self._cache_list))
            df_list = [pd.DataFrame(self._cache_list)]
            if self._cache_df is not None:
                df_list.append(self._cache_df)
                # load cache file again in case it has been updated by other jobs
                try:
                    copyFileEOS(self.cache_fullpath, self.cachefile)
                    df_list.append(pd.read_hdf(self.cachefile, key=self.md5))
                    logger.info('Reloaded cache from %s' % self.cache_fullpath)
                except KeyError:
                    raise
                except Exception:
                    pass
            df = pd.concat(df_list).drop_duplicates(['event', 'jetidx'])
            try:
                df.to_hdf(self.cachefile, key=self.md5, complevel=7, complib='blosc')
                copyFileEOS(self.cachefile, self.cache_fullpath)
                logger.info('New cache file saved to %s' % self.cache_fullpath)
            except Exception:
                logger.error(traceback.format_exc())
            if os.path.exists(self.cachefile):
                os.remove(self.cachefile)

    def predict_with_cache(self, taginfo_producer, event_idx, jet_idx, jet=None):
        outputs = None
        if self._cache_df is not None:
            outputs = self._cache_dict.get((event_idx, jet_idx))
        if outputs is None:
            taginfo = taginfo_producer.load(event_idx)
            outputs = self.predict_one(taginfo, event_idx - taginfo_producer._uproot_start, jet_idx, jet=jet)
            self._cache_list.append({'event': event_idx, 'jetidx': jet_idx, **outputs})
        return outputs


if __name__ == '__main__':
    import time
    import uproot
    import argparse
    parser = argparse.ArgumentParser('TEST')
    parser.add_argument('-i', '--input')
    parser.add_argument('-m', '--model')
    parser.add_argument('-p', '--preprocess')
    parser.add_argument('--make_baseline', action='store_true')
    args = parser.parse_args()

    p = ParticleNetTagInfoMaker()
    tree = uproot.open(args.input)['Events']
    table = tree.arrays(['FatJet*', 'PFCands*', 'SV*'], namedecode='utf-8')
    start = time.time()
    taginfo = p.convert(table)
    diff = time.time() - start
    print('--- Convert inputs: %f s total, %f s per jet ---' % (diff, diff / taginfo['pfcand_mask'].counts.sum()))
#     for k in taginfo:
#         print(k, taginfo[k])

#     jetmass = tree.array('FatJet_msoftdrop')
#     eval_flags = (jetmass > 50) * (jetmass < 200)
#     jetmass = jetmass[eval_flags]
    eval_flags = None

    start = time.time()
    nn = ParticleNetJetTagsProducer(args.model, args.preprocess)
    diff = time.time() - start
    print('--- Setup model: %f s total' % (diff,))

    start = time.time()
    outputs = nn.predict(taginfo, eval_flags)
    diff = time.time() - start
    print('--- Run prediction: %f s total, %f s per jet ---' % (diff, diff / outputs['probQCDbb'].counts.sum()))
#     print(outputs)
#     for k in outputs:
#         print(k, outputs[k].content.mean())

    if 'FatJet_ParticleNetMD_probXbb' in table:
        print('Compare w/ stored values')
        print('Stored values:\n ...', table['FatJet_ParticleNetMD_probXbb'][:5])
        print('Computed values:\n ...', outputs['probXbb'][:5])
        print('Diff (50%, 95%, 99%, 100%) = ', np.percentile(
            np.abs(outputs['probXbb'] - table['FatJet_ParticleNetMD_probXbb']).content, [50, 95, 99, 100]))

#     assert(np.array_equal(jetmass.counts, outputs['probQCDbb'].counts))
    alloutputs = awkward.JaggedArray.zip(outputs)
    if args.make_baseline:
        with open('baseline.awkd', 'wb') as fout:
            awkward.save(fout, alloutputs)
    else:
        if os.path.exists('baseline.awkd'):
            with open('baseline.awkd', 'rb') as fin:
                baseline = awkward.load(fin)
            print("Comparison to baseline:", (alloutputs == baseline).all().all())

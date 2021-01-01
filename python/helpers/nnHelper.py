def convert_prob(jet, sigs, bkgs=None, prefix=''):
    if not jet:
        return -1

    def get(name):
        if isinstance(jet, dict):
            return jet[name]
        else:
            return getattr(jet, name)

    if bkgs is None:
        bkgs = [prefix + n for n in ['QCDbb', 'QCDb', 'QCDcc', 'QCDc', 'QCDothers']]
    else:
        if not isinstance(bkgs, (list, tuple)):
            bkgs = [bkgs]
        bkgs = [prefix + n for n in bkgs]
    bkgsum = sum([get(name) for name in bkgs])

    if sigs is None:
        return bkgsum
    else:
        if not isinstance(sigs, (list, tuple)):
            sigs = [sigs]
        sigs = [prefix + n for n in sigs]
    sigsum = sum([get(name) for name in sigs])

    try:
        return sigsum / (sigsum + bkgsum)
    except ZeroDivisionError:
        return -1


def ensemble(outputs, func):
    keys = outputs[0].keys()
    return {k: func([o[k] for o in outputs]) for k in keys}

def convert_prob(jet, sigs, bkgs=None, prefix=''):
    if not jet:
        return -1

    if bkgs is None:
        bkgs = [prefix + n for n in ['QCDbb', 'QCDb', 'QCDcc', 'QCDc', 'QCDothers']]
    else:
        if not isinstance(bkgs, (list, tuple)):
            bkgs = [bkgs]
        bkgs = [prefix + n for n in bkgs]
    bkgsum = sum([getattr(jet, name) for name in bkgs])

    if sigs is None:
        return bkgsum
    else:
        if not isinstance(sigs, (list, tuple)):
            sigs = [sigs]
        sigs = [prefix + n for n in sigs]
    sigsum = sum([getattr(jet, name) for name in sigs])

    try:
        return sigsum / (sigsum + bkgsum)
    except ZeroDivisionError:
        return -1

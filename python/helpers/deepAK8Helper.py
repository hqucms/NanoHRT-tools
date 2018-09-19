
_nominal_qcd = ('nnQCDbb', 'nnQCDcc', 'nnQCDb', 'nnQCDc', 'nnQCDothers')
_nominal_score_dict = {
    'TvsQCD': (('nnTbcq', 'nnTbqq'), _nominal_qcd),
    'WvsQCD': (('nnWcq', 'nnWqq'), _nominal_qcd),
    'ZvsQCD': (('nnZbb', 'nnZcc', 'nnZqq'), _nominal_qcd),
    'HbbvsQCD': (('nnHbb',), _nominal_qcd),
    }

_decorr_qcd = ('nnMDQCDbb', 'nnMDQCDcc', 'nnMDQCDb', 'nnMDQCDc', 'nnMDQCDothers')
_decorr_score_dict = {
    'TvsQCD': (('nnMDTbcq', 'nnMDTbqq'), _decorr_qcd),
    'WvsQCD': (('nnMDWcq', 'nnMDWqq'), _decorr_qcd),
    'ZHbbvsQCD': (('nnMDZbb', 'nnMDHbb'), _decorr_qcd),
    'bbvsLight': (('nnMDZbb', 'nnMDHbb', 'nnQCDbb'), ('nnQCDcc', 'nnQCDb', 'nnQCDc', 'nnQCDothers')),
    }


def _get_score(jet, score, decorr):
    if not jet:
        return -1
    signames, bkgnames = _decorr_score_dict[score] if decorr else _nominal_score_dict[score]
    sigsum = sum([getattr(jet, name) for name in signames])
    bkgsum = sum([getattr(jet, name) for name in bkgnames])
    try:
        return sigsum / (sigsum + bkgsum)
    except ZeroDivisionError:
        return -1


def get_nominal_score(jet, score):
    return _get_score(jet, score, False)


def get_decorr_score(jet, score):
    return _get_score(jet, score, True)

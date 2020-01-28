def _pass_trig(event, trig):
    # return True only if trig exist and is 1
    try:
        if getattr(event, trig):
            return True
    except:
        pass
    return False

def passTrigger(event, trig_names):
    if not isinstance(trig_names, (list, tuple)):
        trig_names = [trig_names]
    for trig in trig_names:
        if _pass_trig(event, trig):
            return True
    return False



import math

def triangle(t, period=2*math.pi, amplitude=2.0):
    a, p = amplitude, period
    return 4*a/p * abs((((t-p/4)%p)+p)%p - p/2) - a


def sawtooth(t, period=2*math.pi, amplitude = 2):
    return amplitude*(t/period - math.floor(0.5 + t/period))

def square(t, period=2*math.pi, amplitude=2):
    pass

def stairs():
    pass

def ramp():
    pass

def linspace():
    pass

def points(s_ref, steps):
    if isinstance(steps, int):
        steps = [steps]*(len(s_ref) - 1)
    import numpy as np
    diff = np.diff(s_ref,axis=0)
    s = np.array([n/stepi*diff[i]+s_ref[i] for i,stepi in enumerate(steps) for n in range(1,stepi+1)])
    return s

#def cycles(number, fidelity=2, shape="triangle"):
#    return getattr(__module__

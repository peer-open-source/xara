from . import units, mks

def __getattr__(name):
    return units.Dimension(getattr(mks, name))


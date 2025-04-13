from . import units, ips

def __getattr__(name):
    return units.Dimension(getattr(ips, name))


import sys
from . import units
HELP = """
python -m libunit <dst>/<frm>...
                  <sys> <dst>...

Examples

python -m libunit si inch foot pound

python -m libunit m/inch newton/pound
"""
def parse_args(args):
    config = {"a": None, "b": None}
    argi = iter(args[1:])
    for arg in argi:
        if False:
            pass
        if config["a"] is None:
            config["a"] = arg
        elif config["b"] is None:
            config["b"] = arg
    return config

if __name__ == "__main__":
    config = parse_args(sys.argv)
    print(getattr(units.UnitHandler(config["a"]), config["b"]))


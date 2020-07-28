import os
import sys
import re
import datetime as dt


def print_version(modulefile):
    """Print the version of a m.swim module, with the grass version if an active."""
    try:
        import grass.script as grass
        gv = grass.version()
        print("grass v%s, %s" % (gv['version'], gv['build_date']))
    except ImportError:
        pass
    time = dt.datetime.fromtimestamp(os.path.getmtime(__file__))
    with open(modulefile) as f:
        module, version = re.findall(r'm\.swim\.(\w+)\sv(\d+\.\d)', f.read(500))[0]
    print("m.swim.%s v%s, %s" % (module, version, time))
    return sys.exit(0)

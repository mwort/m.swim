"""Input/output related functionality of m.swim.*"""
import numpy as np


def write_csv(filepath, array, column_names, min_width=9, float_columns=[], float_precision=3):
    """A numpy csv writer for float or int arrays."""
    # formate for structure file
    colwidth = [max(i, min_width) for i in map(len, column_names)]
    heacolfmt = ['%%-%ss' % i for i in colwidth] 
    datcolfmt = ['%%%s.%sf' % (i, float_precision) if c in float_columns else '%%%si' % i
                 for i, c in zip(colwidth, column_names)]
    # write out structure file
    with open(filepath, 'w') as strf:
        # header
        strf.write((', '.join(heacolfmt) + '\n') % tuple(column_names))
        # data
        np.savetxt(strf, array, fmt=', '.join(datcolfmt))
    return
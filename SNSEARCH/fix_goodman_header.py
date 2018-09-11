#!/usr/bin/env python

import os
import sys

from ccdproc import CCDData

def fix_goodman_header():

    list_of_files = sys.argv[1:]

    for filename in list_of_files:

        ccd = CCDData.read(filename, unit="adu")

        path, name = os.path.split(filename)
        ## filename = os.path.join(path, 'h' + name)
        ## ccd.write(filename, overwrite=True)

        ccd.write(name, overwrite=True)

if __name__ == '__main__':
    fix_goodman_header()

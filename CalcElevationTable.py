#!/usr/bin/env python3
"""create a table of dB suppression vs elevation"""

from __future__ import print_function
from builtins import range

import sys
import os
import argparse

import pandas as pd
import numpy as np

from calc_width import calc_width

#===== Main ==============================================

def main():
    """
    Calculate the cone width depending on frequency, height and desired dB suppression.
    Height is specified in km above the moon reference radius of 1737.1km.
    Frequency is specified in kHz.
    """

    MHz2KHz = 1000.0
    freq_MHz = 50; #diffraction is worse at lower frequencies, so use low end of band, 50MHz
    freq_kHZ = freq_MHz*MHz2KHz
    height = 0.010; #assume and elevation of 10m (TODO use LRO digital elevation model)

    cone_data = list()
    dBspace = np.linspace(-10, -90, num=160)
    for dB in dBspace:
        #call Neil Bassett's modeling function
        cone_width, cone_err_upper, cone_err_lower = calc_width(freq_kHZ, height, dB)
        row_entry = [freq_MHz, height, dB, cone_width, cone_err_upper, cone_err_lower]
        cone_data.append(row_entry)
        print(row_entry)

    df = pd.DataFrame(cone_data, columns=['Frequency (MHz)', 'Height (km)', 'suppression (dB)', 'cone_width (deg)', "err_upper (deg)", "err_lower (deg)"])
    df.to_csv("./elevation_table.csv")



#===== Official entry point ==============================

if __name__ == '__main__':
    main()
    sys.exit(0)

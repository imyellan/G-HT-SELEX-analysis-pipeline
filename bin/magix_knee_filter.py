#!/usr/bin/env python

import pandas as pd
# import numpy as np
from kneed import KneeLocator
from scipy import stats
import matplotlib.pyplot as plt
import sys

if len(sys.argv) >= 2:
    in_tbl = sys.argv[1]
    skip_header = "infer"

if len(sys.argv) == 3:
    skip_header = None

peaks = pd.read_table(in_tbl, sep='\t', header=skip_header)
coeffs_sorted = peaks.iloc(4).sort_values(ascending=True).values

# remove very low values (less than 1)
coeffs_sorted_thresh = coeffs_sorted[coeffs_sorted >= 1]
ranks = stats.rankdata(coeffs_sorted_thresh, method= 'dense')

kl = KneeLocator(x = ranks, y = coeffs_sorted_thresh, S = 30,
                 curve = "convex", direction = "increasing")
coeff_thresh = kl.knee_y

# filter peaks by coefficient threshold
c17_peaks_filtered = peaks[peaks.iloc(4) >= coeff_thresh]

# save filtered peaks
c17_peaks_filtered.to_csv("knee_filt_peaks.bed", sep='\t', index=False)
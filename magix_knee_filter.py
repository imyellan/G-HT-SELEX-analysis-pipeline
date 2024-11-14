#!/usr/bin/env python

import pandas as pd
# import numpy as np
from kneed import KneeLocator
from scipy import stats
# import matplotlib.pyplot as plt
import sys

in_tbl = sys.argv[1]

peaks = pd.read_table(in_tbl, sep='\t')
coeffs_sorted = peaks["coefficient.br"].sort_values(ascending=True).values
ranks = stats.rankdata(coeffs_sorted, method= 'dense')

kl = KneeLocator(x = ranks, y = coeffs_sorted, S = 10,
                            curve = "convex", direction = "increasing")
coeff_thresh = kl.knee_y

# filter peaks by coefficient threshold
c17_peaks_filtered = peaks[peaks["coefficient.br"] >= coeff_thresh]

# save filtered peaks
c17_peaks_filtered.to_csv("knee_filt_peaks.bed", sep='\t', index=False)
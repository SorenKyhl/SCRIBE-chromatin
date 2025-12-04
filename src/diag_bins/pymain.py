
import numpy as np
import pandas as pd
import bisect

def binDiagonal(d, bins):
    bins = bins
    bin_id =  pd.cut([d], bins, include_lowest=False)#, labels=range(len(bins)-1))
    return bin_id


def generate_lookup(bins, nbeads):
    lookup = []
    curr = 0
    bin_id = 0
    for i in range(nbeads):
        if i >= bins[curr]:
            bin_id += 1
            curr += 1
        lookup.append(bin_id)
    return lookup

if __name__ == '__main__':
    bins = np.loadtxt("bins.txt")
    ans = []
    for i in range(10):
        ans.append(bisect.bisect_right(bins, i))
    print(ans)


    bins = np.loadtxt("lowerbins.tx")
    ans = []
    for i in range(10):
        ans.append(bisect.bisect_left(bins, i))
    print(ans)

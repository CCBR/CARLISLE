#!/usr/bin/env python3
import pandas
import sys
# arguments
# @Inputs
# @param1 $1 = spike-in len file
# @param2 $2 = idx stats file

spikeinlen = pandas.read_csv(sys.argv[1],sep="\t",header=None)
idxstats = pandas.read_csv(sys.argv[2],sep="\t",header=None)

nreads_spikein = sum(list(idxstats[idxstats[0].isin(list(spikeinlen[0]))][2]))

print(nreads_spikein)

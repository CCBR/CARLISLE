#!/usr/bin/env python3
import pandas
import sys
import subprocess
import yaml
# arguments
# @Inputs
# @param1 = R1 fastq file
# @param2 = spike-in len file
# @param3 = genome len file
# @param4 = raw idx stats file
# @param5 = filtered idx stats file - no_dedup
# @param6 = filtered idx stats file - dedup
# @param7 = out yaml file

stats=dict()

r1fastqgz = sys.argv[1]

spikeinlen = pandas.read_csv(sys.argv[2],sep="\t",header=None)
genomelen = pandas.read_csv(sys.argv[3],sep="\t",header=None)
raw_idxstats = pandas.read_csv(sys.argv[4],sep="\t",header=None)
nodedup_idxstats = pandas.read_csv(sys.argv[5],sep="\t",header=None)
dedup_idxstats = pandas.read_csv(sys.argv[6],sep="\t",header=None)

stats["raw_nreads_genome"] = sum(list(raw_idxstats[raw_idxstats[0].isin(list(genomelen[0]))][2]))
stats["raw_nreads_spikein"] = sum(list(raw_idxstats[raw_idxstats[0].isin(list(spikeinlen[0]))][2]))

stats["nodedup_nreads_genome"] = sum(list(nodedup_idxstats[nodedup_idxstats[0].isin(list(genomelen[0]))][2]))
stats["nodedup_nreads_spikein"] = sum(list(nodedup_idxstats[nodedup_idxstats[0].isin(list(spikeinlen[0]))][2]))

stats["dedup_nreads_genome"] = sum(list(dedup_idxstats[dedup_idxstats[0].isin(list(genomelen[0]))][2]))
stats["dedup_nreads_spikein"] = sum(list(dedup_idxstats[dedup_idxstats[0].isin(list(spikeinlen[0]))][2]))

# get number of reads
cmd = "zcat " + r1fastqgz + " | wc -l "
ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
nlines = ps.communicate()[0]
nlines = int(nlines)
npairs = nlines/4
stats["nreads"] = int(npairs*2)

with open(sys.argv[7], 'w') as file:
    dumped = yaml.dump(stats, file)
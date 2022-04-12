#!/usr/bin/env python3

import argparse,pandas,os,yaml,subprocess,uuid
from operator import index
from email import header

randbed=str(uuid.uuid4())+".bed"

parser = argparse.ArgumentParser(description='create bed')
parser.add_argument('--results', required=True, type=str, help="results bed file")
parser.add_argument('--fdr_cutoff', required=True, type=float, help="FDR cutoff")
parser.add_argument('--log2FC_cutoff', required=True, type=float, help="log2FC cutoff")
parser.add_argument('--elbowyaml', required=True, type=str, help="elbow limits")
parser.add_argument('--bed', required=True, type=str, help="output bed with colors")
args = parser.parse_args()

df = pandas.read_csv(args.results,sep="\t",header=0,usecols=['seqnames','start','end','log2FoldChange','padj'])
df['color'] = ["160,160,160"]*len(df.index)
df['name'] = ["."]*len(df.index)
df['strand'] = ["+"]*len(df.index)
df['score'] = [0]*len(df.index)
df['seven'] = [0]*len(df.index)
df['eight'] = [0]*len(df.index)

df.loc[df['log2FoldChange'] < 0,'strand'] = "-"
with open(args.elbowyaml) as f:
    elbowdf = yaml.safe_load(f)
elbowdown = elbowdf['low_limit']
elbowup = elbowdf['up_limit']
if abs(elbowdown) < args.log2FC_cutoff:
    df.loc[(df['padj'] < args.fdr_cutoff) & (df['log2FoldChange'] < elbowdown),'color']="229,255,204"
df.loc[(df['padj'] < args.fdr_cutoff) & (df['log2FoldChange'] < args.log2FC_cutoff),'color']="128,255,0"
if elbowup < args.log2FC_cutoff:
    df.loc[(df['padj'] < args.fdr_cutoff) & (df['log2FoldChange'] > elbowup),'color']="255,153,51"
df.loc[(df['padj'] < args.fdr_cutoff) & (df['log2FoldChange'] > args.log2FC_cutoff),'color']="255,51,51"
df = df.reindex(columns=['seqnames','start','end','name','score','strand','seven','eight','color'])

df.to_csv(randbed,sep="\t",header=False,index=False)

cmd = "sort -S 40G -T /dev/shm -k1,1 -k2,2n "+randbed+" > "+args.bed
s = subprocess.run(cmd,shell=True,check=True,text=True,capture_output=True)

os.remove(randbed)

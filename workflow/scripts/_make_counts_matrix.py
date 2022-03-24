#!/usr/bin/env python3
# This scripts takes inputs.txt as input which has group/bed/bedgraph as tsv columns
# First we concat all bed files and then merge the resultant bed file to get a single bed file
# for regions of interest or ROI
# Next, using the ROI we count the aggregate signal in each ROI for each sample and report the
# sums to countsmatrix.tsv
# sampleinfo.tsv is also created to be used as colData for downstream DESeq2 analysis

import subprocess,argparse,sys,pandas,os,functools,uuid
parser = argparse.ArgumentParser(description='create counts matrix')
parser.add_argument('--bedbedgraph', required=True, type=str, help="list of peak calls and scaled bedgraph file, tab delimited, group|sampleName|bedfile|bedgraph per line")
parser.add_argument('--tmpdir', required=True, type=str, help="TMP dir")
parser.add_argument('--countsmatrix', required=True, type=str, help="output counts matrix TSV")
parser.add_argument('--sampleinfo', required=True, type=str, help="output sample info TSV")

args = parser.parse_args()

randstr = str(uuid.uuid4())

if not os.path.exists(args.tmpdir):
	try:
		os.mkdir(args.tmpdir)
	except:
		exit(args.tmpdir + ": Folder doesnt exist and cannot be created!")

if not os.access(args.tmpdir,os.W_OK):
	exit(args.tmpdir + ": Folder is not writeable!")

bedbedgraph = pandas.read_csv(args.bedbedgraph,header=None,sep="\t")
bedbedgraph.columns=["group","sampleName","bed","bedgraph"]

concatfile = os.path.join(args.tmpdir,randstr+".concat.bed")
cmd = "cat " + " ".join(list(bedbedgraph["bed"])) + " | cut -f1-3 | sort -k1,1V -k2,2n > " + concatfile
concat = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

mergedfile = os.path.join(args.tmpdir,randstr+".merged.bed")
cmd = "bedtools merge -i " + concatfile + " > " + mergedfile
merge = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

sampleinfofile = open(args.sampleinfo,'w')
sampleinfofile.write("samplename\tgroup\n")
counts=dict()
for i,row in bedbedgraph.iterrows():
	group = row['group']
	bg = row['bedgraph']
	samplename = row['sampleName']
	filepath,filename = os.path.split(bg)
	basename,ext = os.path.splitext(filename)
	sampleinfofile.write("%s\t%s\n"%(samplename,group))
	cmd = "bedtools intersect -wa -wb -a "+bg+" -b "+mergedfile
	bt = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
	peakCounts=dict()
	for l in bt.stdout.split("\n"):
		l=l.strip().split("\t")
		if len(l) != 7:
			continue
		peakID = l[4] + ":" + l[5] + "-" + l[6]
		score = (int(l[2]) - int(l[1])) * float(l[3])
		if not peakID in peakCounts:
			peakCounts[peakID]=0.0
		peakCounts[peakID]+=score
	counts[samplename] = pandas.DataFrame.from_dict(peakCounts,orient='index',columns=[samplename])
sampleinfofile.close()
counts_matrix_df = functools.reduce(lambda x,y:x.merge(y,left_index=True,right_index=True),counts.values())
counts_matrix_df.to_csv(args.countsmatrix,sep="\t",header=True,index=True,index_label="peakID")
#os.rmdir(args.tmpdir)
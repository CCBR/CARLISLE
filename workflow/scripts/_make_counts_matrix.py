#!/usr/bin/env python3
# This scripts takes inputs.txt as input which has group/bed/bedgraph as tsv columns
# First we concat all bed files and then merge the resultant bed file to get a single bed file
# for regions of interest or ROI
# Next, using the ROI we count the aggregate signal in each ROI for each sample and report the
# sums to countsmatrix.tsv
# sampleinfo.tsv is also created to be used as colData for downstream DESeq2 analysis

import subprocess,argparse,sys,pandas,os,functools,uuid
parser = argparse.ArgumentParser(description='create counts matrix')
parser.add_argument('--bedbedgraph', required=True, type=str, help="list of peak calls and scaled bedgraph file, tab delimited, group|sampleName|bedfile|bedgraph|scalingFactor|fragmentsBed per line")
parser.add_argument('--tmpdir', required=True, type=str, help="TMP dir")
parser.add_argument('--countsmatrix', required=True, type=str, help="bedgraph-AUC-based output counts matrix TSV")
parser.add_argument('--fragmentscountsmatrix', required=True, type=str, help="fragmentsBed-based output counts matrix TSV")
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
bedbedgraph.columns=["group","sampleName","bed","bedgraph","scalingFactor","fragmentsBed"]

concatfile = os.path.join(args.tmpdir,randstr+".concat.bed")
cmd = "cat " + " ".join(list(bedbedgraph["bed"])) + " | cut -f1-3 | sort -k1,1 -k2,2n > " + concatfile
concat = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

mergedfile = os.path.join(args.tmpdir,randstr+".merged.bed")
cmd = "bedtools merge -i " + concatfile + " > " + mergedfile
merge = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

sampleinfofile = open(args.sampleinfo,'w')
sampleinfofile.write("samplename\tgroup\n")
counts=dict() # bedgraph AUC based counts
fcounts=dict() # fragment based counts
for i,row in bedbedgraph.iterrows():
	group = row['group']
	bg = row['bedgraph']
	samplename = row['sampleName']
	sf = row['scalingFactor']
	fbed = row['fragmentsBed']
	filepath,filename = os.path.split(bg)
	basename,ext = os.path.splitext(filename)
	sampleinfofile.write("%s\t%s\n"%(samplename,group))
	cmd = "bedtools intersect -wa -wb -a "+bg+" -b "+mergedfile
	print("Now running: %s"%(cmd))
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
	cmd = "bedmap --echo-ref-name --count "+mergedfile+" "+fbed
	print("Now running %s"%(cmd))
	bo = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
	fPeakCounts=dict()
	for l in bo.stdout.split('\n'):
		l = l.strip().split('|')
		if len(l)!=2: 
			continue
		fPeakCounts[l[0]]=l[1]
	fcounts[samplename] = pandas.DataFrame.from_dict(fPeakCounts,orient='index',columns=[samplename])
sampleinfofile.close()
counts_matrix_df = functools.reduce(lambda x,y:x.merge(y,left_index=True,right_index=True),counts.values())
counts_matrix_df.to_csv(args.countsmatrix,sep="\t",header=True,index=True,index_label="peakID")
fcounts_matrix_df = functools.reduce(lambda x,y:x.merge(y,left_index=True,right_index=True),fcounts.values())
fcounts_matrix_df.to_csv(args.fragmentscountsmatrix,sep="\t",header=True,index=True,index_label="peakID")

#os.rmdir(args.tmpdir)
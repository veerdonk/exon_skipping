#! /usr/bin/python3

import math
import sys

outDir = "jobs/eb/{}/".format(sys.argv[1])
baseName = sys.argv[1]
nJobs = math.ceil(2014/100)
pathToDir = "/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/"


nStudies = 2011
ranges = list()
for i in range(0, nStudies, 100):
	if(i+100 <= nStudies):
		ranges.append([i+1, i+100])
	else:
		ranges.append([i+1, nStudies])

	

for i in range(nJobs):
	outfile = open(outDir + baseName + str(i) + ".sh", "w")
	jobName = str(i) + baseName

	outfile.write('''#!/bin/bash
#SBATCH --job-name={}
#SBATCH --output={}.out
#SBATCH --error={}.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=35gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load RPlus/3.4.4-foss-2015b-v18.03.2
export R_LIBS="/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/R_packages/:$R_LIBS"

Rscript {}find_other_gene_cmd.R \\
{}{} \\
{}{} \\
{} \\
{} {} {}'''.format(jobName, jobName, jobName, \
pathToDir, pathToDir+"output/eb/{}/".format(sys.argv[1]) ,jobName + "out.csv", pathToDir, \
"bed/EB_exons.bed", ranges[i][0], ranges[i][1],\
sys.argv[2], sys.argv[3]))


	outfile.close()



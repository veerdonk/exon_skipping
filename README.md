# Analyzing naturally ocurring exon skipping in COL7A1

### About the project
The code for my graduation project "Analyzing naturally ocurring exon skipping in COL7A1" is stored in this repository. During this project I investigated the possibility of using online RNA-seq data repositories like recount2 to find skipped exons in COL7A1. The data I collected was then used to generate predictive data on the ease of artificially skipping exons using AONs.

If you are interested in the results or the methods used please contact me so I can provide you with a copy of the report.

### Project structure
The scripts used for the exon skipping analysis can be found in the cluster_scripts folder. When using this script the path to the folder you want the studies to be saved will have to be edited. In the home directory of this repository two other scripts are available. One can be used to download all the neccessary data from recount2 and the other is used to generate figures for the data produced by the exon skipping analysis.

### Prerequisites
##### Python 3.6.x
Python is used to run createJobScripts.py and will have to installed. A download of the most recent version can be found [here](https://www.python.org/downloads/)
##### R
The scripts used for the exon skipping analysis are written in R and require R as well as several packages to be installed. R itself can be downloaded [here](https://www.r-project.org/).
The packages used are:
* recount
* ggplot2
* gridExtra
* plyr
* data.table
These will have to be installed prior to running the scripts as they will not function without these packages.

##### Computing cluster
While not a strict necessity the use of a computing cluster can drastically bring down the time needed to process all the data in recount2. createJobScripts.py generates job files that can be used on a slurm cluster like UMCG boxy.

##### Data
Before any analysis can be run some data needs to be collected. For your gene of interest you need:
* Gene name (symbol)
* Ensembl ID
* UCSC ID
* .bed file containing exon annotation
* (optional) downloaded data from recount

### Usage
Firstly createJobScripts.py has to be used to generate the shell scripts that start each individual job. Before running the script the location of the .bed file containing the exon annotation of all/some genes has to be added to the script. The syntax for running this script is: 

`python3 createJobScripts.py [gene_name] [Ensembl_ID] [UCSC_ID]`

An example for COL7A1 would look like this: 

`python3 createJobScripts.py COL7A1 ENSG00000114270.16 uc003ctz.3`

NOTE: the Ensembl ID version has to match the ID version stored in recount.
This script takes one or two seconds to run and produces ~20 shell scripts in a directory named after the gene.

Once these batch files have been produced the jobs can be started by running a shell command like:

`for i in $( ls *.sh ); do sbatch --qos=regular $i; done`

This starts a job for each (~20) .sh file in the current directory with the urgency set to regular. The shell scripts will then call the main exon skipping script (find_other_gene_cmd.R) which in turn will start the analysis. If all data from recount2 has been downloaded prior to the run the total runtime per job should be around 45-60 minutes, this is recommended as you will not be dependant on an internet connection and the servers of recount2. The output will be ~20 .csv files all containing the exons as columns and the read counts of the samples of 100 studies as rows

To download all of the data necessary to analyze exon skipping using recount2 data download_all_studies.R can be used. Inside the script the destination folder has to be changed to your preferred directory. running this script from the commandline will then start downloading all of the studies available on recount2 (2004 at time of writing) along with the .bed file containing junction start and stop locations. Running this script will obviously require some disk space but will also speed up any future runs significantly.

To produce some figures plot_all.R can be used. This script cannot be used from the commandline. The easiest way to use it is to open it in RStudio and run the part of the script that produce the figures you want/need.

### Descriptions
createJobScripts.py - Takes some identifiers and produces a number of shell scripts that can be used to start jobs on a slurm cluster.
download_all_studies.R - Downloads all necessary data from recount2 so they will not have to be downloaded each run.
find_other_gene_cmd.R - The main exon skip finding script. This is called upon by the shell scripts created by createJobScripts.py and produces ~20 .csv output files.
plot_all.R - Script that takes the output of find_other_gene_cmd.R, parses it and then uses the data to produce a number of insightfull figures.

### Author
David van de Veerdonk

### Acknowledgements
Peter van den Akker - Instructor UMCG & project assignment
Freerk van Dijk - Instructor UMCG
Joeri van der Velde - Instructor UMCG
Dave Langers - Intsructor Hanze

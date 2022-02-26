
# RNApipe 
## 1. Installation
### 1.1 Install Miniconda
As this workflow is based on the workflow management system snakemake and conda.
We strongly recommend installing miniconda3 with python3.
For specific installation methods and usage methods of RNApipe, please refer to the RNApipe_Documentation.pdf.
###  1.2 Install RNApipe
Clone the repository:
```
git clone https://github.com/ywu019/RNApipe.git
```
Create the environment:
```
conda env create -n RNApipe -f envs/envs.yaml
```
Activate the environment:
```
conda activate RNApipe
```

## 2. Prepare input files
Several input files are required in order to run the workflow, a genome file (.fa), an annotation file (.gff/.gtf) and compressed sequencing files (.fastq.gz).
| File type | Description |
| ------ | ------ |
| genome.fa | user-provided genome file containing the genome sequence |
| annotation.gtf | user-provided annotation file with genomic features |
| <condition>_<replicate>_<end>.fastq.gz | user-provided compressed sequencing files |
| config.yaml | configuration file to customize the workflow |
| samples.tsv | sample file describing the relation between the input fastq files |
### 2.1  Annotation.gff and genome.fa
We recommend retrieving both the genome and the annotation files for your organism from National Center for Biotechnology Information (NCBI)
Note: if you use custom annotation files, ensure that you adhere to the gtf standard
 
### 2.2 Input .fastq files
These are the input files provided by the user. Both single end and paired end data is supported.
Note: Please ensure that you compress your files in .gz format and <condition>_<replicate>_<end>.fastq.gz (root_1_R1.fastq.gz if paired end data /root_1.fastq.gz if single end data )

### 2.3 Set up configuration
Modify the metafile describing your data configs/samples.tsv .
| condition | replicate |
| ------ | ------ |
| root | 1 |
| root | 2 |
| root | 3 |
| shoot | 1 |
| shoot | 2 |
| shoot | 3 |

Customize the workflow based on your need in configs/config.yaml .It contains the following variables:

* PROJECT: Project name
* READSPATH: The path to fastq files
* SAMPLES: configs/samples.tsv
* END: sequencing paired-end or single-end
* OUTPUTPATH: Thw path for final outputs
* GENOME: The path of genome files
* ANNOTATION: The path of annotation files
* CONTROL: Control group in comparison
* TREAT: Processing groups in comparison
* SPEICE_ABBREVIATION: Abbreviation of species name,Query in the species_abbreviation.tsv

## 3. Run RNApipe
```
python main.py
```

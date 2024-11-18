[![DOI](https://zenodo.org/badge/829822589.svg)](https://zenodo.org/doi/10.5281/zenodo.13347495)

# NanoMnT
A collection of Python scripts for (1) correcting sequencing errors within STR regions from Oxford Nanopore (ONT) sequencing data and (2) genotyping monoallelic STR regions. 

## Installation

Before installing, ensure that the following packages are installed:
  
  - Matplotlib >= 3.7.1
  - Numpy >= 1.20.3
  - Pysam >= 0.20.0
  - Pandas >= 2.0.0
  - Scipy >= 1.7.1
  - Seaborn >= 0.13.0

After installing these packages, use the following commands to install NanoMnT. 
*Feel free to create a new conda environment for NanoMnT.
```
git clone https://github.com/18parkky/NanoMnT.git
cd NanoMnT
pip install .
```

NanoMnT requires [Krait](https://github.com/lmdu/krait) output files as inputs, which needs to be installed separately.  
Krait output files for the human genome (T2T-CHM13v2.0) is included in NanoMnT:```NanoMnT/ref/krait.t2tchm13.tsv```.
NanoMnT has been tested in Python version 3.9.7. and should work in later versions.

## Commands & Parameters

NanoMnT has three main commands: **getAlleleTable**, **getLocusTable**, and **findInformativeLoci** as well as several utility commands.

### getAlleleTable
The ```getAlleleTable``` command requires an input BAM and Krait file (to specify STR loci of interest) as inputs, and produces a single TSV file, namely, the allele table. The allele table contains information about every reads that aligned to any of the given STR locus.
```
usage: getAlleleTable.py [-h] -b PATH_BAM -s PATH_STR_TSV -r PATH_REF_GENOME [-m MAPPING_QUALITY] [-f FLANKING] [-rf REALIGNMENT_FLANKING] [-t THREADS] [--sc] [-out DIR_OUT]
```
**Required parameters**:
- **-b**: PATH to input BAM
- **-s**: PATH to Krait output  
- **-r**: PATH to reference genome
   
**Optional parameters**:
- -m: Minimum MAPQ when processing reads. Reads with MAPQ below this will be discarded (default: 60)
- -f: Length of flanking sequences to be written on allele table (default: 12)
- -rf: Length of flanking sequences of STR to use as reference, during re-alignment process (default: 20000)
- -t: Number of threads to use for multiprocessing (default: 4)
- --sc: Option that tells NanoMnT to look for CB and UMI for each read (default: False)
- -out: Directory to write output files (default: current directory)

### getLocusTable
getLocusTable takes the generated allele table as input and summarizes information for each detected STR locus, generating a single TSV file called the locus table.  
```
usage: getLocusTable.py [-h] -at PATH_STR_ALLELE_TABLE -read_selection READ_SELECTION [--get_allele_size_info] [-cov MIN_COVERAGE] [-t THREADS] [-out DIR_OUT]
```
**Required parameters**:
- **-at**: PATH to allele table
- **-read_selection**: method to use when selecting which set of reads to use when estimating STR allele size:
    - <ins>all_reads</ins>: use all reads
    - <ins>fw_rv_reads</ins>: use forward reads or reverse reads, depending on the repeat unit of the STR. Currently, this applies only to A-/T-/AC-/TG-repeats. All STR with other repeat units will be subjected to use all reads.
    - <ins>ML</ins>: similar to fw_rv_reads, but determine whether to use forward reads or reverse reads via ML prediction of sequencing accuracy using STR flanking sequences. Read NanoMnT paper for more context.

      
**Optional parameters**:
- --get_allele_size_dist_info: option that tells NanoMnT to generate additional files (default: False):  
  (1) a PNG file for allele histogram (red dashed line represents the genotyped allele, blue dashed line represents the reference allele. If genotyped allele == reference allele, only one line will be visible, as shown below. )
  
    ![alt](https://github.com/18parkky/NanoMnT/blob/main/example/Ax12_chr1_MATERNAL_19288456-19288467.png)  

  (2) allele frequency table   
    | Allele | Frequency    |  
    | :-----: | :---: |  
    | 9 | 0.1   |  
    | 10 | 0.1   |  
    | 11 | 0.15   |
    | 12 | 0.5   |  
    | 13 | 0.15   |
  
- -ref: directory of the reference genome. This parameter is required when specifying ML mode as read_selection parameter (default: None)
- -cov: minimum required coverage when processing a given STR locus. Loci with coverage lower than this value will not be processed (default: 20)
- -t: number of threads to use for multiprocessing (default: 4)
- -out: Directory to write output file to (default: current directory)

### findInformativeLoci
getLocusTable takes the generated allele table as input and summarizes information for each detected STR locus, generating a single TSV file called the locus table.  
findInformativeLoci takes allele table and locus table of both tumor and paired normal sample as inputs, and searches for informative STR loci. i.e., STR loci that differ in tumor allele and normal allele. High number of informative STR loci may indicate MSI-H phenotype.
```
usage: findInformativeLoci.py [-h] -nat PATH_NORMAL_AT -nlt PATH_NORMAL_LT -tat PATH_TUMOR_AT -tlt PATH_TUMOR_LT [-cov MIN_COVERAGE] [-n_plots NUM_PLOTS] [-out DIR_OUT]
```
**Required parameters**:
- **-nat**: PATH to allele table of the normal sample
- **-nlt**: PATH to locus table of the normal sample
- **-tat**: PATH to allele table of the tumor sample
- **-tlt**: PATH to locus table of the tumor sample

**Optional parameters**:
- -cov: Minimum coverage required when analyzing STR locus. Loci with coverage lower than this value will not be processed (default: 20)
- -n_plots: Generate allele histogram of top __ informative loci (default: 10)
- -out: Directory to write output files (default: current directory)


### utility script: filterBAMbyMAPQ
NanoMnT also has a utility command: **filterBAMbyMAPQ**, which takes a BAM file as input and filters out non-Primary reads with user-specifiable MAPQ threshold. Although not necessary, you can first preprocess your BAM file with this command before running getAlleleTable & getLocusTable.  
```
usage: filterBAMbyMAPQ.py [-h] -b PATH_INPUT_BAM -m MAPPING_QUALITY [-P] [--remove_secondary_reads] [--remove_supplementary_reads] [-out DIR_OUT]
```
**Required parameters**:
- **-b**: PATH to input BAM
- **-m**: MAPQ threshold to filter reads. Note that the maximum MAPQ output by minimap2 is 60
**Optional parameters**:
- -P: Only output primary reads
- --remove_secondary_reads: Remove secondary reads
- --remove_supplementary_reads: Remove supplementary reads
- -out: Directory to write output file to (default: current directory)
  
## Output files

By default, running each command will each generate one file (+one log file).  

### getAlleleTable → allele table
As mentioned above, getAlleleTable produces the allele table, of which each row contains the following information:
- read name
- locus
- repeat unit
- allele
- STR size of reference
- STR size of read
- base quality of STR and adjacent (12 nt for each direction) bases
- left flanking sequence of STR
- right flanking sequence of STR
- Cell barcode 
- Unique molecular identifier
- Corrected allele (error-corrected allele)
- Flag of the read

### getLocusTable → locus table
Similarly, a successful run of getLocusTable also generates a single tab-separated file, referred to as the locus table. 
Locus table contains information about each STR region that have been detected within your input BAM.
Each row of locus table contains the following information:
- STR locus (e.g., chr1:10001-100012)
- Major allele<sup>1</sup>	
- Repeat unit (e.g., CAG)
- STR size of reference
- Effective number of allele in the the given STR locus<sup>2</sup>	
- Percentage of reference allele in the given STR locus
- Correction rate<sup>3</sup>
- Predicted number of allele in the given STR locus<sup>4</sup>	
- **Genotyped allele**
- Peak prominence of the most frequently observed allele<sup>5</sup>
  
1) <ins>Major allele</ins> is defined by the most frequently found allele within the allele size histogram of your data.
2) <ins>Effective number of allele</ins> is not meant to be interpreted directly; rather, it describes the dispersion of the allele size histogram. If the STR locus is believed to be monoallelic, low effective number of allele would be a good indication.
3) <ins>Correction rate</ins> represents the percentage of reads that have been corrected during ```getAlleleTable```.
4) <ins>Predicted number of allele</ins> is measured by counting the number of prominent peaks in the allele size histogram of the given STR locus.
5) <ins>Peak prominence</ins> is calculated using SciPy's find_peaks function. Furthermore, the most frequent allele does not necessarily equal the genotyped allele (ONT can generate misleading STR allele size histogram, and the most frequently observed allele is often not the correct allele). Nevertheless, this metric often provides good QC metric, as <ins>high peak prominence most often indicates good STR size estimation</ins>. 

### findInformativeLoci → informative loci table
A successful run of findInformativeLoci generates a single tab-separated file, referred to as informative loci table.
Each row of informative loci table contains the following information:
- STR locus
- Repeat unit
- STR size of reference
- Coverage of normal sample
- Coverage of tumor sample
- Peak prominence of normal STR allele size histogram
- Peak prominence of tumor STR allele size histogram
- Histogram similarity: similarity between STR allele size histogram between normal vs. tumor sample
- Locus score, calculated by: (1-histogram similarity)*peak prominence of normal sample.

Among these, the most important information is the locus score. High locus score indicates highly informative STR locus. 

### filterBAMbyMAPQ → Filtered BAM  
filterBAMbyMAPQ generates a filtered, sorted BAM file along with the BAM index file.  

## Tutorial - Identifying MSI status of cancer cell lines using SG-NEx dataset
Although ONT sequencing suffers heavily in STR regions, it can still be used to detect MSI status in cancer sequencing data, nonetheless.
We provide a tutorial (~/tutorial/tutorial.ipynb) for detecting MSI status by inspecting STR allele size profile of cancer ONT data.  

  
## Citation  
If you use NanoMnT in your research, please cite this paper.  


NanoMnT is under active development, so feel free to make suggestions! 

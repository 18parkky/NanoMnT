[![DOI](https://zenodo.org/badge/829822589.svg)](https://zenodo.org/doi/10.5281/zenodo.13347495)

# Table of contents
1. [Installation](#installation)<br/>
2. [Commands & Parameters](#commands--parameters)<br/>
3. [Output files](#output-files)<br/>
4. [Tutorial - Identifying MSI status of cancer cell lines using SG-NEx dataset](#tutorial---identifying-msi-status-of-cancer-cell-lines-using-sg-nex-dataset)<br/>
5. [Citation](#citation)<br/>

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
Krait output files for 3 versions the human genome (GRCh37, GRCh38, and T2T-CHM13v2.0) are included in NanoMnT GitHub repository under:`~/NanoMnT/ref/krait/`.
Krait configuration for minimum repeats for each STR type is as follows: mono-: 10, di-: 7, tri-: 5, tetra-: 4, penta-: 4, and hexa-: 4. 
NanoMnT has been tested in Python version 3.9.7. and should work in later versions.

## Commands & Parameters
NanoMnT has four main commands: `estimateLociCoverage`, `getAlleleTable`, `getLocusTable`, and `findInformativeLoci` as well as several utility commands.

### `estimateLociCoverage`
Note: The `estimateLociCoverage` command is actually closer to a utility command rather than the main command, as it's job is to estimate coverage of each locus and filter out STR loci with insufficient coverage, **prior to running `getAlleleTable`**. Nonetheless, unless (1) your sequencing data evenly covers the entire genome (e.g., whole-genome sequencing) or (2) you have thoughtfully selected your STR loci by manually pre-processing Krait outputs, it is recommended that you run this to save time because `getAlleleTable` can take some time.

Given a BAM file and a Krait output file, this script estimates the coverage of each locus and annotates the Krait file with this information. Two Krait files will be saved to disk: (1)
`~/{filename}.coverageLabeled.tsv` and (2) `~/{filename}.coveragedLabeled.highCov.tsv`. Using the latter Krait file for getAlleleTable will significantly reduce runtime.
```
usage: estimateLociCoverage.py [-h] -b PATH_BAM -s PATH_STR_TSV [-cov MIN_COVERAGE] [-out DIR_OUT]
```
**Required parameters**:
- `-b`: PATH to input BAM
- `-s`: PATH to Krait output  
   
**Optional parameters**:
- `-cov`: Minimum coverage required when analyzing STR locus. Loci with coverage lower than this value will not be processed (default: 20)
- `-out`: Directory to write output files (default: current directory)

### `getAlleleTable`
The `getAlleleTable` command requires an input BAM and Krait file (to specify STR loci of interest) as inputs, and produces a single TSV file, namely, the allele table. The allele table contains information about every reads that aligned to any of the given STR locus.
```
usage: getAlleleTable.py [-h] -b PATH_BAM -s PATH_STR_TSV -r PATH_REF_GENOME [-m MAPPING_QUALITY] [-f FLANKING] [-rf REALIGNMENT_FLANKING] [-t THREADS] [--sc] [-out DIR_OUT]
```
**Required parameters**:
- `-b`: PATH to input BAM
- `-s`: PATH to Krait output  
- `-r`: PATH to reference genome
   
**Optional parameters**:
- `-m`: Minimum MAPQ when processing reads. Reads with MAPQ below this will be discarded (default: 60)
- `-f`: Length of flanking sequences to be written on allele table (default: 12)
- `-rf`: Length of flanking sequences of STR to use as reference, during re-alignment process (default: 20000)
- `-t`: Number of threads to use for multiprocessing (default: 4)
- `--sc`: Option that tells NanoMnT to look for CB and UMI for each read (default: False)
- `-out`: Directory to write output files (default: current directory)

### `getLocusTable`
`getLocusTable` takes the generated allele table as input and summarizes information for each detected STR locus, generating a single TSV file called the locus table.  
```
usage: getLocusTable.py [-h] -at PATH_STR_ALLELE_TABLE -read_selection READ_SELECTION [--get_allele_size_info] [-cov MIN_COVERAGE] [-t THREADS] [-out DIR_OUT]
```
**Required parameters**:
- `-at`: PATH to allele table
- `-read_selection`: method to use when selecting which set of reads to use when estimating STR allele size:
    - <ins>all_reads</ins>: use all reads
    - <ins>fw_rv_reads</ins>: use forward reads or reverse reads, depending on the repeat unit of the STR. Currently, this applies only to A-/T-/AC-/TG-repeats. All STR with other repeat units will be subjected to use all reads.
    - <ins>ML</ins>: similar to fw_rv_reads, but determine whether to use forward reads or reverse reads via ML prediction of sequencing accuracy using STR flanking sequences. Read NanoMnT paper for more context.

      
**Optional parameters**:
- --get_allele_size_dist_info: option that tells NanoMnT to generate additional files (default: False):  
  (1) a PNG file for allele histogram (red dashed line represents the genotyped allele, blue dashed line represents the reference allele. If `genotyped allele == reference allele`, only one line will be visible, as shown below. )
  
    ![alt](https://github.com/18parkky/NanoMnT/blob/main/example/Ax12_chr1_MATERNAL_19288456-19288467.png)  

  (2) allele frequency table   
    | Allele | Frequency    |  
    | :-----: | :---: |  
    | 9 | 0.1   |  
    | 10 | 0.1   |  
    | 11 | 0.15   |
    | 12 | 0.5   |  
    | 13 | 0.15   |
  
- `-ref`: directory of the reference genome. This parameter is required when specifying ML mode as read_selection parameter (default: None)
- `-cov`: minimum required coverage when processing a given STR locus. Loci with coverage lower than this value will not be processed (default: 20)
- `-t`: number of threads to use for multiprocessing (default: 4)
- `-out`: Directory to write output file to (default: current directory)

### `findInformativeLoci`
`findInformativeLoci` takes allele table and locus table of both tumor and paired normal sample as inputs, and searches for informative STR loci. i.e., STR loci that differ in tumor allele and normal allele. High number of informative STR loci may indicate MSI-H phenotype.
```
usage: findInformativeLoci.py [-h] -nat PATH_NORMAL_AT -nlt PATH_NORMAL_LT -tat PATH_TUMOR_AT -tlt PATH_TUMOR_LT [-cov MIN_COVERAGE] [-n_plots NUM_PLOTS] [-out DIR_OUT]
```
**Required parameters**:
- `-nat`: PATH to allele table of the normal sample
- `-nlt`: PATH to locus table of the normal sample
- `-tat`: PATH to allele table of the tumor sample
- `-tlt`: PATH to locus table of the tumor sample

**Optional parameters**:
- `-cov`: Minimum coverage required when analyzing STR locus. Loci with coverage lower than this value will not be processed (default: 20)
- `-n_plots`: Generate allele histogram of top __ informative loci (default: 10)
- `-out`: Directory to write output files (default: current directory)

### `filterBAMbyMAPQ` (utility script)
NanoMnT also has a utility script: `python ~/NanoMnT/NanoMnT/utility/filterBAMbyMAPQ.py` (please note that this is a Python script, not a command that can be directly used from the terminal), which takes a BAM file as input and filters out non-Primary reads with user-specifiable MAPQ threshold. Although not necessary, you can first preprocess your BAM file with this command before running `getAlleleTable` & `getLocusTable`.  
```
usage: filterBAMbyMAPQ.py [-h] -b PATH_INPUT_BAM -m MAPPING_QUALITY [-P] [--remove_secondary_reads] [--remove_supplementary_reads] [-out DIR_OUT]
```
**Required parameters**:
- `-b`: PATH to input BAM
- `-m`: MAPQ threshold to filter reads. Note that the maximum MAPQ output by `minimap2` is 60
**Optional parameters**:
- `-P`: Only output primary reads
- `--remove_secondary_reads`: Remove secondary reads
- `--remove_supplementary_reads`: Remove supplementary reads
- `-out`: Directory to write output file to (default: current directory)

## Output files

**We recommend going through the tutorial notebook(s). It will help users understand the outputs & how to interpret them**
By default, running each command will each generate one file (+one log file).  

### `estimateLociCoverage` → Two Krait files
`estimateLociCoverage` produces two Krait files: (1) the input Krait file labeld with estimated coverage without any filtering and (2) a smaller Krait file where STR loci with insufficient coverage (<`cov`) are filtered out.
Uusing the second Krait file for `getAlleleTable` will significantly reduce runtime.

### `getAlleleTable` → allele table
As mentioned above, `getAlleleTable` produces the 'allele table', of which each row contains the following information:
- Read name (`read_name`)
- Locus (`locus`)
- Repeat unit (`repeat_unit`)
- STR allele sequence reported by the read (e.g., AAAAAAAA) (`allele`)
- STR allele (=number of repeats) of the reference genome (`reference_STR_allele`)
- Left flanking sequence of STR (`left_flanking_seq`)
- Right flanking sequence of STR (`right_flanking_seq`)
- Flag of the read (`flag`)
- Cell barcode (`CB`)
- Unique molecular identifier (`UMI`)
- Corrected allele (error-corrected allele) (`corrected_allele`)
- Editing distance between the uncorrected allele and the corrected allele (`editing distance`)
- STR allele reported by the read (`read_STR_allele`)

### `getLocusTable` → locus table
Similarly, a successful run of `getLocusTable` also generates a single tab-separated file, referred to as the 'locus table'. 
Locus table contains information about each STR region that have been detected within your input BAM.
Each row of locus table contains the following information:
- Chromosome (`chromosome`)
- Start position of the STR locus in the reference genome (`start`)
- End position of the STR locus in the referenc genome (`end`)
- STR locus name (e.g., chr1:10001-100012) (`locus`)
- Repeat unit (`repeat_unit`)
- STR allele (=number of repeats) of the reference genome (`reference_STR_allele`)
- Percentage of the reads that report the reference allele (`percentage_of_reference_allele`)
- Percentage of reads that have been corrected (`correction_rate`)
- Number of reads that have been corrected (`corrected_read_count`)
- <ins>**STR allele estimated by NanoMnT**</ins> (`allele`)
- Position of the most prominent peak in the allele histogram (to visualize the allele histogram, use `--get_allele_size_info`) (`peak_pos`)
- Prominence of the peak<sup>1</sup> (`peak_prominence`)
- Number of peaks detected by `scipy.signal.find_peaks` (`num_peaks`) (unless the sample is MSI, high `num_peaks` indicates poor sequencing accuracy or high ploidy)
- STR allele relative to the reference genome (`int_relative_allele_size`) (e.g., v>0 indicates that the allele is longer than the reference allele)
- Area of allele histogram relative to the reference (`histogram_area_relative_to_reference`)<sup>2</sup>
- Read selection method used (the `-read_selection` parameter) (`read_selection`)

1) <ins>Peak prominence</ins> is calculated using SciPy's find_peaks function. Furthermore, the most frequent allele does not necessarily equal the genotyped allele (ONT can sometimes generate misleading STR allele size histogram, and the most frequently observed allele is often not the correct allele). Nevertheless, this metric often provides good QC metric, as <ins>high peak prominence most often indicates good STR size estimation</ins>.<br/>
2) This metric is calculated as the sum of `(STR allele - reference allele) × (observed frequency of the allele)` for each observed allele. Very similar to `int_relative_allele_size` in terms of interpretation.
   
### `findInformativeLoci` → informative loci table
A successful run of `findInformativeLoci` generates a single tab-separated file, referred to as informative loci table.
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

### `filterBAMbyMAPQ` → Filtered BAM  
filterBAMbyMAPQ generates a filtered, sorted BAM file along with the BAM index file.  

## Tutorial - Identifying MSI status of cancer cell lines using SG-NEx dataset
Although ONT sequencing suffers heavily in STR regions, it can still be used to detect MSI status in cancer sequencing data, nonetheless.
We provide a tutorial (~/tutorial/tutorial.ipynb) for detecting MSI status by inspecting STR allele size profile of cancer ONT data.  

  
## Citation  
If you use NanoMnT in your research, please cite this paper.  


NanoMnT is under active development, so feel free to make suggestions! 

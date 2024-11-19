import os
import time
import pysam
import logging
import argparse
import numpy as np
import pandas as pd

import NanoMnT.nanomnt_utility as nanomnt_utility

def main():
    start_time  = time.time()

    script_description = "Given a BAM file and a Krait output file, this script estimates the coverage of each locus and annotates the Krait file with this information. Two Krait files will be saved to disk: (1) ~/{filename}.coverageLabeled.tsv and (2) ~/{filename}.coveragedLabeled.highCov.tsv. Using the latter Krait file for getAlleleTable will significantly reduce runtime."
    parser = argparse.ArgumentParser(description=script_description)

    # Read paramters
    # Required parameters
    parser.add_argument('-b', '--PATH_bam',
                        help="PATH to input BAM file (must be sorted and indexed)",
                        required=True)
    parser.add_argument('-s', '--PATH_str_tsv',
                        help="PATH to STR list file generated using either Krait or Pytrf (.tsv)",
                        required=True,
                        )
    parser.add_argument('-cov', '--min_coverage',
                        help="Minimum coverage required for each locus. Loci with coverage below this value will be excluded from output Krait file (default: 20)",
                        required=False, type=int, default=20)
    parser.add_argument('-out', '--DIR_out',
                        help='Directory to write output files (default: current directory)',
                        required=False,
                        type=str,
                        default=os.getcwd(),
                        )
    
    args = vars(parser.parse_args())

    PATH_bam        = args["PATH_bam"]
    PATH_str_tsv    = args["PATH_str_tsv"]
    cov_threhold    = args['min_coverage']
    DIR_out         = args["DIR_out"]
    
    bam_filename = os.path.splitext(os.path.basename(PATH_bam))[0]
    krait_filename = os.path.splitext(os.path.basename(PATH_str_tsv))[0]

    # Create log file
    nanomnt_utility.checkAndCreate(DIR_out)
    PATH_log = f"{DIR_out}/nanomnt.estimateLociCoverage.{bam_filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')
        
    # Load files
    krait = pd.read_csv(PATH_str_tsv, sep='\t')
    bamfile = pysam.AlignmentFile(PATH_bam, 'rb')
    
    # Estimate coverage
    logging.info(f'Estimating coverage of {len(krait)} loci')
    col_estimated_cov = list()
    for tup in krait.itertuples():
        chrom, start, end = tup.sequence, tup.start, tup.end
        middle_pos = int( (start + end) / 2 )
        estimated_cov = 0 
        for pileupcolumn in bamfile.pileup( chrom, middle_pos, middle_pos+1, max_depth=10000, min_base_quality=0, truncate=True ):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_refskip: estimated_cov += 1
        col_estimated_cov.append(estimated_cov)
    krait['estimated_coverage'] = col_estimated_cov
    logging.info(f'Finished estimation\t(elapsed time: {nanomnt_utility.getElapsedTime(start_time)} seconds)')
    
    krait_covered = krait[krait['estimated_coverage']!=0].copy()
    num_covered_loci = len(krait_covered)
    logging.info(f'{num_covered_loci} loci ({round(100*num_covered_loci/len(krait), 2)}%) were covered by input BAM')
    logging.info(f'Mean and std of estimated coverage among these loci:\t{np.mean(krait_covered["estimated_coverage"])}, {np.std(krait_covered["estimated_coverage"])}')

    krait_high_covered = krait_covered[krait_covered['estimated_coverage']>=cov_threhold].copy()
    logging.info(f'Number of loci with coverage ≥ {cov_threhold}: {len(krait_high_covered)}')
    
    krait_high_covered.reset_index(inplace=True, drop=True)

    # Save to disk
    krait.to_csv(f'{DIR_out}/{krait_filename}-{bam_filename}.coverageLabeled.tsv', sep='\t', index=False)
    krait_high_covered.to_csv(f'{DIR_out}/{krait_filename}-{bam_filename}.coverageLabeled.highCov.tsv', sep='\t', index=False)

    logging.info(f"Finished estimateLociCoverage.py\t(Total time taken: {nanomnt_utility.getElapsedTime(start_time)} seconds)")

if __name__ == "__main__":
    main()
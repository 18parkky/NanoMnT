import gc
import os
import time
import math 
import glob
import pysam
import logging
import argparse
import subprocess
import numpy as np
import Levenshtein
import pandas as pd
import multiprocessing

import seaborn as sns 
import matplotlib.pyplot as plt

import NanoMnT.nanomnt_utility as nanomnt_utility

def runFindInformativeLoci( PATH_normal_at, PATH_normal_lt, 
                           PATH_tumor_at, PATH_tumor_lt, 
                           min_coverage, num_plots, 
                           start_time, PATH_log, output_filename, DIR_out ):
    
    logging.basicConfig(filename=PATH_log, level=logging.INFO)

    df_normal_at    = pd.read_csv(PATH_normal_at, sep='\t')
    df_normal_lt    = pd.read_csv(PATH_normal_lt, sep='\t')
    df_tumor_at     = pd.read_csv(PATH_tumor_at, sep='\t')
    df_tumor_lt     = pd.read_csv(PATH_tumor_lt, sep='\t')
    
    # (1) Filter out loci with insufficient coverage

    set_sufficiently_covered_normal_loci    = set()
    set_sufficiently_covered_tumor_loci     = set()

    for tup in df_normal_lt.itertuples():
        if tup.corrected_read_count >= min_coverage:
            set_sufficiently_covered_normal_loci.add( tup.locus )
            
    for tup in df_tumor_lt.itertuples():
        if tup.corrected_read_count >= min_coverage:
            set_sufficiently_covered_tumor_loci.add( tup.locus )

    set_intersecting_loci = set_sufficiently_covered_normal_loci.intersection( set_sufficiently_covered_tumor_loci ) 

    df_normal_at    = df_normal_at[(df_normal_at["locus"].isin( set_intersecting_loci ))]
    df_tumor_at     = df_tumor_at[(df_tumor_at["locus"].isin( set_intersecting_loci ))]
    
    
    # (2) Calculate histogram distance between Normal vs. Tumor
    dict_locus_to_alleleHistograms = dict() 
    read_selection = df_tumor_lt.iloc[0].read_selection
    logging.info(f"Read selection used for locus table: {read_selection}")
    if read_selection not in ["all_reads", "fw_rv_reads"]:
        logging.error(f"Invalid read selection method detected for locus table! Try running getLocusTable.py again!\t(elapsed time: {nanomnt_utility.getElapsedTime(start_time)} seconds)")


    for locus, edf in df_normal_at.groupby("locus"):
        if read_selection == "all_reads":
            edf_oi = edf 
        elif read_selection == "fw_rv_reads":
            edf_oi = edf[( edf["flag"]==nanomnt_utility.dict_repeatUnit_to_favorableFlag[ edf.iloc[0].repeat_unit ] )]
            
        dict_locus_to_alleleHistograms[locus] = [nanomnt_utility.getAlleleHistogram( edf_oi[["repeat_unit", "reference_STR_allele", "corrected_allele"]].dropna(), nanomnt_utility.allele_histogram_width )]

    for locus, edf in df_tumor_at.groupby("locus"):
        if read_selection == "all_reads":
            edf_oi = edf 
        elif read_selection == "fw_rv_reads":
            edf_oi = edf[( edf["flag"]==nanomnt_utility.dict_repeatUnit_to_favorableFlag[ edf.iloc[0].repeat_unit ] )]
            
        dict_locus_to_alleleHistograms[locus].append( nanomnt_utility.getAlleleHistogram( edf_oi[["repeat_unit", "reference_STR_allele", "corrected_allele"]].dropna(), nanomnt_utility.allele_histogram_width ) )
        
    dict_locus_to_histogramSimilarity = dict()
    for locus, alleleHistograms in dict_locus_to_alleleHistograms.items():
        dict_locus_to_histogramSimilarity[locus] = 1 - nanomnt_utility.calcHistogramDistance( alleleHistograms[0], alleleHistograms[1], "cosine" )
        
    # (3) Get other locus info from LocusTables (e.g., repeat unit, reference allele, peak prominence, etc.)
    df_normal_lt = df_normal_lt[(df_normal_lt["locus"].isin( dict_locus_to_histogramSimilarity.keys() ))]
    df_tumor_lt  = df_tumor_lt[(df_tumor_lt["locus"].isin( dict_locus_to_histogramSimilarity.keys() ))]
    
    dict_locus_to_LocusInfo = dict()

    for tup in df_normal_lt.itertuples():
        dict_locus_to_LocusInfo[tup.locus] = [ tup.repeat_unit, tup.reference_STR_allele, tup.corrected_read_count, tup.peak_prominence ]
    for tup in df_tumor_lt.itertuples():
        dict_locus_to_LocusInfo[tup.locus].insert( 3, tup.corrected_read_count )
        dict_locus_to_LocusInfo[tup.locus].append( tup.peak_prominence )
        dict_locus_to_LocusInfo[tup.locus].append( dict_locus_to_histogramSimilarity[tup.locus] )
    
    df_informativeLoci_info = pd.DataFrame.from_dict(dict_locus_to_LocusInfo, orient='index', columns=['repeat_unit', 'reference_STR_allele', 'normal_cov', 'tumor_cov', 
                                                                                                       'normal_peak_prominence', 'tumor_peak_prominence', 'histogram_similarity'])
    
    # (4) Calculate locus score (score represents how much you can trust the result)
    list_loci_score = list()

    for tup in df_informativeLoci_info.itertuples():
        score = ( 1 - tup.histogram_similarity ) * ( tup.normal_peak_prominence )
        list_loci_score.append( score )

    df_informativeLoci_info["locus_score"] = list_loci_score
    df_informativeLoci_info.sort_values('locus_score', inplace=True, ascending=False)
    df_informativeLoci_info.reset_index(inplace=True, names="locus")
    
    logging.info(f"{len(df_informativeLoci_info)} loci were both covered in normal and tumor samples")
    
    # (5) Generate STR allele size histogram plots
    if num_plots > 0:
        logging.info(f"Creating STR allele size histogram plots for top informative {num_plots} loci")
        DIR_histogram_out = f"{DIR_out}/{output_filename}_histograms"
        nanomnt_utility.checkAndCreate( DIR_histogram_out )

        plot_count = 0
        for tup in df_informativeLoci_info.dropna().itertuples():
        
            f = sns.lineplot( dict_locus_to_alleleHistograms[tup.locus][0], alpha=0.5, color="seagreen" )
            f = sns.scatterplot( dict_locus_to_alleleHistograms[tup.locus][0], alpha=0.5, color="seagreen" )
            
            f = sns.lineplot( dict_locus_to_alleleHistograms[tup.locus][1], alpha=0.5, color="darkred" )
            f = sns.scatterplot( dict_locus_to_alleleHistograms[tup.locus][1], alpha=0.5, color="darkred" )
            
            f.set(title=f"{tup.locus}\n{tup.repeat_unit}x{tup.reference_STR_allele}_locusScore:{round(tup.locus_score, 1)}")
            f.set_xlabel("STR allele")
            f.set_ylabel("Percentage (%)")

            fig = f.get_figure()
            fig.savefig( f"{DIR_histogram_out}/{plot_count+1}_{tup.locus}.png" )
            plt.clf()

            plot_count += 1
            if plot_count == num_plots:
                break 
        
    # (6) Write to disk and finish
    df_informativeLoci_info.to_csv(f"{DIR_out}/{output_filename}.informativeLoci.tsv", sep='\t', index=False)
    return 

def main():
    start_time  = time.time()

    github_link = "https://github.com/18parkky/NanoMnT"
    script_description = f"Given paired allele table & locus table for paired normal and tumor samples, search for informative loci. For more information, see GitHub: {github_link}"
    parser = argparse.ArgumentParser(description=script_description)

    # Required parameters
    parser.add_argument('-nat', '--PATH_normal_at',      
                        help="PATH to paired normal allele table", 
                        required=True,
                        )
    
    parser.add_argument('-nlt', '--PATH_normal_lt',      
                        help="PATH to paired normal locus table", 
                        required=True,
                        )
    
    parser.add_argument('-tat', '--PATH_tumor_at',      
                        help="PATH to tumor allele table", 
                        required=True,
                        )
    
    parser.add_argument('-tlt', '--PATH_tumor_lt',      
                        help="PATH to tumor locus table", 
                        required=True,
                        )
    
    # Optional parameters
    parser.add_argument('-cov', '--min_coverage',       
                        help="Minimum coverage required for analyzing STR locus. Loci with coverage lower than this value will not be processed (default: 20)", 
                        required=False, 
                        type=int, 
                        default=20
                        )
    
    parser.add_argument('-n_plots', '--num_plots',       
                        help="Generate allele histogram of top __ informative loci (default: 10)", 
                        required=False, 
                        type=int, 
                        default=10
                        )
    
    parser.add_argument('-out', '--DIR_out',          
                        help='Directory to write output files (default: current directory)', 
                        required=False, 
                        type=str, 
                        default=os.getcwd() 
                        )
    
    args = vars(parser.parse_args())

    PATH_normal_at   = args["PATH_normal_at"]
    PATH_normal_lt   = args["PATH_normal_lt"]
    PATH_tumor_at    = args["PATH_tumor_at"]
    PATH_tumor_lt    = args["PATH_tumor_lt"]
    min_coverage    = args["min_coverage"]
    num_plots       = args["num_plots"]
    DIR_out         = args["DIR_out"]
    
    normal_filename = os.path.splitext( os.path.basename( PATH_normal_at ))[0]
    tumor_filename  = os.path.splitext( os.path.basename( PATH_tumor_at ))[0]
    output_filename = f"{normal_filename}_{tumor_filename}"

    nanomnt_utility.checkAndCreate( DIR_out )
    
    PATH_log = f"{DIR_out}/nanomnt.findInformativeLoci.{output_filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing parameters:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')

    runFindInformativeLoci( PATH_normal_at, PATH_normal_lt, 
                           PATH_tumor_at, PATH_tumor_lt, 
                           min_coverage, num_plots, 
                           start_time, PATH_log, output_filename, DIR_out )

    logging.info(f"Finished findInformativeLoci.py\t(Total time taken: {nanomnt_utility.getElapsedTime(start_time)} seconds)")

if __name__ == "__main__":
    main()
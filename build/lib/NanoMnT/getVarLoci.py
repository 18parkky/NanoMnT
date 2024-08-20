import math
import logging
import argparse
import pandas as pd 

import seaborn as sns 
import matplotlib.pyplot as plt 
sns.set_style("darkgrid", {"grid.color": ".6", "grid.linestyle": ":"})

from . import nanomnt_utility

def main():
    start_time  = time.time()

    script_description = f"By inputting allele table and locus table of normal and tumor samples, getVarLoci finds STR loci that have different genotyped results and summarizes them into a TSV file. For more information, see GitHub documentation"
    parser = argparse.ArgumentParser(description=script_description)
    
    script_PATH     = os.path.dirname(os.path.dirname(__file__))

    # Required parameters
    parser.add_argument('-n_at', '--dir_normal_STR_allele_table',   help="Directory of normal STR allele table", required=True)
    parser.add_argument('-n_lt', '--dir_normal_STR_locus_table',    help="Directory of normal STR locus table", required=True)
    parser.add_argument('-t_at', '--dir_tumor_STR_allele_table',    help="Directory of tumor STR allele table", required=True)
    parser.add_argument('-t_lt', '--dir_tumor_STR_locus_table',     help="Directory of tumor STR locus table", required=True)

    # Optional parameters
    parser.add_argument('--get_allele_size_dist_info',  help="Generate allele size distribution information of each locus (default: False)", action='store_true')    
    parser.add_argument('-cov', '--min_coverage',       help="Minimum coverage required for analyzing STR locus (default: 20)", required=False, type=int, default=20)
    parser.add_argument('-PATH', '--PATH_out',          help='PATH of the output files (default: current directory)', required=False, type=str, default=os.getcwd() )

    args = vars(parser.parse_args())
    
    dir_STR_allele_table        = args["dir_STR_allele_table"]
    read_selection              = args["read_selection"]
    dir_ref_genome              = args["dir_ref_genome"]
    get_allele_size_dist_info   = args["get_allele_size_dist_info"]
    min_coverage                = args["min_coverage"]
    threads                     = args["threads"]
    PATH_out                    = args["PATH_out"]
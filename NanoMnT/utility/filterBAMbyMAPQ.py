import os 
import time
import pysam
import argparse

def main():
    start_time = time.time()

    script_description = ''

    parser = argparse.ArgumentParser(description=script_description)
    # Required parameters
    parser.add_argument('-b', '--dir_input_bam', help='Directory of input BAM file', required=True)
    parser.add_argument('-m', '--mapping_quality', help='The highest MAPQ score the aligner assigns (minimap2: 60, bowtie2: 42)', required=True, type=int)
    # Optional parameters
    parser.add_argument('-P', '--primary_only', help='Reads that aren\'t primary reads are removed', action='store_true')
    parser.add_argument('--remove_secondary_reads', help='Remove secondary reads', action='store_true')
    parser.add_argument('--remove_supplementary_reads', help='Remove supplementary reads', action='store_true')
    parser.add_argument('-PATH', '--PATH_out',          help='PATH of the output files (default: current directory)', required=False, type=str, default=os.getcwd() )

    args = vars(parser.parse_args())

    dir_bam_file = args['dir_input_bam']
    file_name = os.path.splitext(os.path.basename(dir_bam_file))[0]
    mapq_max = args['mapping_quality']
    only_primary = args['primary_only']
    remove_secondary = args['remove_secondary_reads']
    remove_supplementary = args['remove_supplementary_reads']
    PATH_out    = args["PATH_out"]

    dir_out_bam = f'{PATH_out}/{file_name}.MAPQ_filtered.bam'
    bamfile = pysam.AlignmentFile(dir_bam_file)
    new_bamfile = pysam.AlignmentFile(dir_out_bam, "wb", template=bamfile)

    if only_primary == True:
        for read in bamfile.fetch():
            if (read.mapq >= mapq_max and read.is_secondary == False and read.is_supplementary == False):
                new_bamfile.write(read)
    else:
        if remove_secondary == True and remove_supplementary == True: # Same as using -P option
            for read in bamfile.fetch():
                if (read.mapq >= mapq_max and read.is_secondary == False and read.is_supplementary == False):
                    new_bamfile.write(read) 
        
        # Keep supplementary alignments, remove secondary alignments
        elif remove_secondary == True and remove_supplementary == False:
            for read in bamfile.fetch():
                if (read.mapq >= mapq_max and read.is_secondary == False):
                    new_bamfile.write(read)

        # Keep secondary alignments, remove supplementary alignments
        elif remove_secondary == False and remove_supplementary == True:
            for read in bamfile.fetch():
                if (read.mapq >= mapq_max and read.is_supplementary == False):
                    new_bamfile.write(read) 

        # Keep both secondary and supplementary alignments
        else:
            for read in bamfile.fetch():
                if (read.mapq >= mapq_max):
                    new_bamfile.write(read)
    new_bamfile.close()
    pysam.index(dir_out_bam)

    print(f"Finished filterBAMbyMAPQ.py\t(Total time taken: {round(time.time() - start_time, 2)} seconds)")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
This script removes all variants in homopolymers longer or equal to threshold -n, and the 
variant with a QUAL below -q minscore. 
The homopolymer length retained is the longest one between the reference and the variant.
'''

#--------------------------------------------------------------------------------------------------
#                                        Import modules
#--------------------------------------------------------------------------------------------------
#TODO check all these modules are required

from io import StringIO
import os
from os import path
import argparse
import Bio
from Bio import SeqIO
from .util import wf_parser


#--------------------------------------------------------------------------------------------------
#                                   Command line examples
#--------------------------------------------------------------------------------------------------

'''
/home/aline/scripts/filter_small_variant_file.py \
  -f /home/aline/raw_data/references/GCA_900235035.2_KP_7435-4_genomic_mito_pla.fasta \
  -v /home/aline/analysis/ABL_Lonza/SA2185_AQM424/variant_calling_small/medaka.sorted_test.vcf \
  -t fasta -q 10 -n 10
'''

    #--------------------------------------------------------------------------------------------------
    #                                          functions
    #--------------------------------------------------------------------------------------------------

def get_homopolymer_len(ref_seq, alt_seq, line_split, record_dict):
    ''' return the length of the longest homopolymer (variant or ref)'''
    
    common_part = ''
    continue_homopolymer = True
    keep_variant = True
    
    if len(ref_seq) > len(alt_seq):	

        longest_seq = ref_seq
        shortest_seq = alt_seq
        # print('case 1')
        
    else  :
        longest_seq = alt_seq
        shortest_seq = ref_seq
        # print('case 2')

    # determine the common part
    for index in range(0,len(shortest_seq)):
        if longest_seq[index] == shortest_seq[index]:
            common_part += alt_seq[index]
        else :
            break
                        
    # print(f"common part : {common_part}")
            
    # we take the first different letter as a potential part of an homopolymer
    homopoly_letter = longest_seq[len(common_part)]
    # print("homopolymer letter : ",homopoly_letter)
    
    homopoly_len = 1
    # we check if we can increase the homopolymer in the different letters between ref and alt
    for index in range((len(common_part)+1),(len(longest_seq))):
        if longest_seq[index] == homopoly_letter : 
            homopoly_len +=1
        else : 
            continue_homopolymer = False
            break
        # print("index : ", index, "letter : ",longest_seq[index], "length : ",homopoly_len)
                        
    print("homopolymer len in variant file : ",homopoly_len)
    # now we need to check if the homopolymer continues in the sequences
    if continue_homopolymer :
        seqname = line_split[0]
        
        # position of the first ref base
        position = int(line_split[1])
                
        # we select the sequence from the ending position of the variant until +25bp after the startin position
        following_sequence = record_dict[seqname].seq[position + len(ref_seq) -1 : position +30]                
        # print("following sequence : ",following_sequence)
                
        # now we see how long is the homopolymer : 
        for base in following_sequence.upper() : 
            if base == homopoly_letter :
                homopoly_len += 1
            else : 
                break     

        # print("final homopolymer length : ", homopoly_len)                        

    return homopoly_len


#--------------------------------------------------------------------------------------------------
#                                                main
#--------------------------------------------------------------------------------------------------

def main(args):

    minscore = int(args.minscore)
    fasta_file_name = args.fasta
    file_type = args.type

    vcf_file_name = args.vcf

    homopolymer_threshold = int(args.homopolymer_threshold)

    # if no output filename is given, build one
    if args.output_file != None :
        output_file_name = args.output_file
    else :
        output_file_name = os.path.splitext(vcf_file_name)[0]+ '_filtered_homopolymers_and_qscore_'+ str(minscore) + '.vcf'


    vcf_file = open(vcf_file_name,'r')
    output_file = open(output_file_name,'w')


    nb_in_homopoly = 0
    nb_qual_too_low = 0
    tot_variants = 0

    record_dict = SeqIO.index(fasta_file_name,file_type)

    # pass the vcf header
    vcf_line = vcf_file.readline()
    if vcf_line:
        while vcf_line[0] == '#' : 
            output_file.write(vcf_line)
            vcf_line = vcf_file.readline()
            
        # read the vcf 
        while vcf_line :
            line_split = vcf_line.split(' ')

            tot_variants += 1
            # print(f"*****\n{vcf_line} ****\n")
            qual = line_split[5]
            # continue only if the quality score of the variant is equal or higher as minscore
            if float(qual) >= minscore : 
                # get the reference and alternative
                ref_seq = line_split[3]
                alt_seq = line_split[4]
                
                keep_variant = True
                
                # continue only if the two are shorter than 6bp
                if len(ref_seq) < 6 and len(alt_seq) < 6 and len(ref_seq) != len(alt_seq) :

                    homopoly_len = get_homopolymer_len(ref_seq, alt_seq, line_split, record_dict)
                    
                    # print("final homopolymer length : ", homopoly_len)                        
                    # we compare the length of the homopolymer with the threshold
                    if homopoly_len >= homopolymer_threshold :
                        keep_variant = False
                        nb_in_homopoly += 1

                    print("keep variant : " , keep_variant)       
                # if we keep the variant we write it
                if keep_variant : 
                    output_file.write(vcf_line)
            
            else :
                nb_qual_too_low += 1
                
            vcf_line = vcf_file.readline()
            
        print("\n*** Result summary ***")
        print(f"{tot_variants} were considered in total")
        print(f"{nb_in_homopoly} variants were removed because they were in homopolymers longer or equal to {homopolymer_threshold}")
        print(f"{nb_qual_too_low} variants were removed because their quality score is lower than {minscore}")


    vcf_file.close()
    output_file.close()

#--------------------------------------------------------------------------------------------------
#                                   Command line parsing
#--------------------------------------------------------------------------------------------------

def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")

    parser.add_argument('-q', '--minscore', help = 'minimum score for a variant to be retained', required = True)
    parser.add_argument('-f','--fasta', help = 'fasta file with the sequences', required = True)
    parser.add_argument('-t','--type', help = 'fasta or fastq', required = True)
    parser.add_argument('-v','--vcf', help = 'variant file in vcf format', required = True)
    parser.add_argument('-n', '--homopolymer_threshold', help = 'threshold after which a variant is filtered out', required = True)
    parser.add_argument('-o','--output_file', help = 'optional, output file name, if not given, the output filename will constructed from the input file name, same folder.', required = False)

    return parser


if __name__ == "__main__":
    args = argparse().parse_args()
    main(args)

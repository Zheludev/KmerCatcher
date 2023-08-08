## this code was originally written by Andy Fire and then modified by Ivan Zheludev
## code to quickly identify reads or other sequences that have bits of illumina linker
    ## works by identifying k-mer matches derived from a database file against a given input multi-fasta

## INZ edits (quality of life) - 07 Aug 2023:
    ## accept multi-line fasta
    ## use in terminal as a script
    ## automatically filter (keep) zero caught sequences

import argparse
from tqdm import tqdm

def main(k_len, database_file, sequence_file, output_file, debug):
    
    ## set variables
    
    #k_len = 12
    #database_file = 'test_1.fa'
    #sequence_file = 'test_2.fa'
    #output_file = 'out.fa'
    #debug = 1
    
    ## parse if output needed
    
    if output_file:
        keep_zero = 1
    else:
        keep_zero = 0
    ##
    
    ## open files and convert to single line
    
    def opener(filename):
        open_file = open(filename, mode='r')
        open_list = open_file.read().splitlines()
        for entry_ind, entry in enumerate(open_list):
            open_list[entry_ind] = entry
        ##
        return open_list
    ##
    
    database_list = opener(database_file)
    sequence_list = opener(sequence_file)
    
    def singleline(in_fasta_list):
        out_fasta_list = []
        new_seq_line = []
        for line in in_fasta_list:
            if line.startswith('>'):
                if new_seq_line:
                    out_fasta_list.append(''.join(new_seq_line))
                    new_seq_line = []
                out_fasta_list.append(line)
            else:
                new_seq_line.append(line.strip())
            ##
        if new_seq_line:
            out_fasta_list.append(''.join(new_seq_line))
        ##
        return out_fasta_list
    ##
    
    database_list = singleline(database_list)
    sequence_list = singleline(sequence_list)
    
    def antisense(s):
        return s.replace('G','c').replace('C','g').replace('A','t').replace('T','a')[::-1].upper()
    ##
    
    ## need to build a dictionary of the database's k-mers to iterate over
        ## need a dummy placeholder value of '0'
    
    database_dict = {}
    for line in database_list:
        if not line.startswith('>'):
            anti_line = antisense(line)
            for index in range(len(line) - k_len + 1):
                database_dict[line[index:index+k_len]] = 1
                database_dict[anti_line[index:index+k_len]] = 1
            ##
        ##
    ##
    
    ## count k-mers that match between database and input
        ## and optionally keep the zero caught input sequences
    
    def counter(in_fasta_element):
        caught = 0
        for index in range(len(in_fasta_element) - k_len + 1):
            in_fasta_kmer = in_fasta_element[index:index + k_len]
            if in_fasta_kmer in database_dict:
                caught += 1
            ##
        ##
        return caught
    ##
    
    if debug:
        print("-----=====-----")
        print("searching for k-mer matches between " + database_file + " (database) and " + sequence_file + " (query)")
        print('seqID', 'k-mers caught', sep='\t')
    ##
    
    keep_list = []
    
    for line_ind, line in tqdm(enumerate(sequence_list)):
        if line.startswith('>'):
            seqID = line
            seq = sequence_list[line_ind + 1]
            count = counter(seq)
            if debug:
                print(seqID, str(int(count)), sep='\t')
            ##
            if keep_zero:
                if (count == 0):
                    keep_list.append(seqID)
                    keep_list.append(seq)
                ##
            ##
            count = 0
        ##
    ##
    
    if debug:
        input_seq_count = str(int((len(sequence_list)/2)))
        output_seq_count = str(int((len(keep_list)/2)))
        print("-----=====-----")
        print("started with " + input_seq_count + " sequences")
        print("identified " + output_seq_count + " sequences with zero k-mer matches")
    ##
    
    def saver(input_name, input_list):
        name_obj = open(input_name, "w")
        for element in input_list:
            if isinstance(element, list):
                element = "\t".join(element)
            name_obj.write(element + "\n")
        ##
        name_obj.close()
    ##
    
    if keep_zero:
        if not keep_list:
            if debug:
                print("no zero k-mer sequences to keep")
            ##
        else:
            if debug:
                print("saving zero k-mer sequences to " + output_file)
            ##
            saver(output_file, keep_list)
        ##
    ##
    
##

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-k_len', type=int, default = 12, help='k-mer length (nt). default = 12')
    parser.add_argument('-db', type=str, help='Database. a .fasta format input file from which k-mers are extracted to match against')
    parser.add_argument('-i', type=str, help='Query. a .fasta format input file in which sequences are compared to the database to find k-mer matches')
    parser.add_argument('-o', type=str, default = '', help='OPTIONAL: Output. a .fasta format output file of all query sequences with zero database matches')
    parser.add_argument('-debug', type=int, default = 1, help='print results. default = 1')
    
    args = parser.parse_args()
    main(args.k_len, args.db, args.i, args.o, args.debug)
##
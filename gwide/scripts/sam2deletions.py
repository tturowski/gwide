#!/usr/bin/env python
import sys, re
import pandas as pd

filepath = sys.argv[1]
df1_temp = pd.read_csv(filepath, sep="\t", names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
df2_comments = df1_temp[df1_temp['QNAME'].str.contains('^@')] #header line start with @
df3_deletions = df1_temp[~df1_temp['QNAME'].str.contains('^@') & df1_temp['CIGAR'].str.contains('D')] #selecting only reads with deletions
df4_single_dels = df3_deletions[df3_deletions['CIGAR'].str.contains('1D')] #only single deletions; will accept 11D and 21D - to be improoved

def del_position(CIGAR = str()):
    '''Takes reference and data DataFrame's

    CIGAR : str
      CIGAR string from SAM file to be parsed
    Returns
    -------
    position of 1st deletion in SAM read
    '''
    CIGAR_elements = re.findall(r'\d{1,2}[A-Z]',CIGAR)
    first_left_match = re.search(r'\d{1,2}M',CIGAR).group(0)
    first_del = re.search(r'\d{1,2}D',CIGAR).group(0)
    #selecting elements between fist MATCH and first DELETION
    range_of_elem = [no for no, elem in enumerate(CIGAR_elements) if elem == first_left_match or elem == first_del]
    #adding
    to_add_elem = CIGAR_elements[range_of_elem[0]:range_of_elem[1]] #taking positions of numbers to pick up range
    nt_to_add = [int(re.search(r'\d{1,2}',i).group(0)) for i in to_add_elem] #extracting numbers
    return sum(nt_to_add)

s1_to_first_del = df4_single_dels['CIGAR'].apply(del_position)

#preparing output file
df5_dels = df4_single_dels.copy()
df5_dels['POS'] = df5_dels['POS'].astype(int) + s1_to_first_del
df5_dels['CIGAR'] = "1D" #mandatory - samtools consider reads without CIGAR as unmapped
df5_dels['SEQ'] = "*"
df5_dels['FLAG'] = df5_dels['FLAG'].astype(int)
df5_dels['MAPQ'] = df5_dels['MAPQ'].astype(int)
df5_dels['PNEXT'] = df5_dels['PNEXT'].astype(int)
df5_dels['TLEN'] = df5_dels['TLEN'].astype(int)

#saving header
df2_comments.to_csv(filepath.replace(".sam","_DEL_ONLY.sam"), index=None, header=None, sep="\t")
#saving reads
with open(filepath.replace(".sam","_DEL_ONLY.sam"), 'a') as output_file:
    df5_dels.to_csv(output_file, index=None, header=None, sep="\t")
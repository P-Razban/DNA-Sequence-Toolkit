from Bio import AlignIO
from Bio import SeqIO
from Bio import pairwise2
import copy
import pandas as pd 

def process_data(csv_path, align_path, output_file):
# Collecting Data
    miRDP_output = pd.read_csv(csv_path)
    align = AlignIO.read(align_path, "clustal")

    total = {}
    
#iterating over each row, extracing and formatting data
    for _, row in miRDP_output.iterrows():
        seeds_location = [int(itm) for itm in row['seeds_location'].strip('][').split(', ')]
        seed_len = [len(row['seeds'].strip('][').split(', ')[i].replace("'", "")) for i in range(len(row['seeds'].strip('][').split(', ')))]
        seq_first_part = row['sequence'].split('.')[0]
        full_seq = row['sequence']
        mr = row['mirna']
        
#finding wild seq and calculating length
        for record in align:
            if seq_first_part in record.id:
                temp_record = record
                temp_len = len(temp_record.seq) - len(temp_record.seq.lstrip('-'))  #first nucleotide position in wilde seq
                
#iterating over other seq and collecting data
                for record in align:
                    seq_A, seq_B, score = [], [], []
                    if not seq_first_part in record.id:
                        total.setdefault('mirna', []).append(mr)
                        total.setdefault('sequence', []).append(full_seq)
                        total.setdefault('seeds_location', []).append(seeds_location)
                        total.setdefault('second_id', []).append(record.id)
                        
#finding seed locations and sequences in varients
                        for i in range(len(seeds_location)):
                            temp_range_record = temp_record.seq[temp_len:temp_len+seeds_location[i]].count('-') #number of gaps within wild seq
                            temp_range = temp_range_record + seeds_location[i] + temp_len - 1                   #first nucleotide of seed 

                            seq1 = temp_record.seq[temp_range: temp_range + seed_len[i]]  #wild seq
                            seq_A.append(seq1)
                            seq2 = record.seq[temp_range: temp_range + seed_len[i]]       #varients seq
                            seq_B.append(seq2)
                            
#pairwise alignment between wild and varient seeds and give them score
                            alignments = pairwise2.align.globalxx(seq1, seq2)
                            score.append(alignments[0].score)

                        total.setdefault('seqA', []).append(seq_A)
                        total.setdefault('seqB', []).append(seq_B)
                        total.setdefault('score', []).append(score)

#reading data and group it by certain columns
    data_frame = pd.DataFrame(total)
    data_frame.to_csv('total_score.csv')

    data_frame = pd.read_csv('total_score.csv')
    sliced_df = data_frame.copy()

    grouped_df = sliced_df.groupby(['sequence', 'mirna'])
    
#finding rows within each group with different 'score' values and collecting their indices
    diff_indices = []
    for _, group in grouped_df:
        if len(group) > 1:
            diff_rows = group['score'].apply(lambda x: x != group['score'].max())
            diff_indices.extend(group[diff_rows].index.tolist())

    diff_indices = list(set(diff_indices))
    final_form_df = sliced_df.loc[diff_indices].sort_index()
    final_form_df.to_csv(output_file, index=False)


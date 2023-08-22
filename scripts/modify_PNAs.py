import sys
from Bio import SeqIO
from difflib import SequenceMatcher
from Bio.SeqIO import SeqRecord
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import pandas as pd

seq_path = sys.argv[1]
res_path = sys.argv[2]

seqs = []
output_df = pd.DataFrame(columns=["ASO", "ASO_seq", "target_seq", "location",
                                      "SC_bases", "pur_perc", "long_pur_stretch", "Tm",
                                      "OT_tot", "OT_TIR"])

for record in SeqIO.parse(seq_path, "fasta"):
    aso = record.seq

    aso_target = record.seq.reverse_complement()
    maxcomp = SequenceMatcher(None, aso, aso_target).find_longest_match(0, len(aso), 0, len(aso_target)).size
    
    aso_name = record.id
    # check for longest purine stretch and purine perc:
    pur = 0
    longest_purine_stretch = 0
    curstretch = 0
    for base in aso.__str__():
        if base in ["A", "G"]:
            curstretch += 1
            pur += 1
            if curstretch > longest_purine_stretch:
                longest_purine_stretch += 1
        else:
            curstretch = 0
    pur_perc = "{:.2f}".format((pur / len(aso.__str__())) * 100)

    added_row = pd.Series([aso_name, aso.__str__(), aso_target.__str__(), None,
                           maxcomp, pur_perc, longest_purine_stretch, None, None, None], index=output_df.columns)
    # output_df = output_df.append(added_row, ignore_index=True)
    output_df = pd.concat([output_df, pd.DataFrame([added_row])], ignore_index=True)
    seqs.append(SeqIO.SeqRecord(aso_target, record.id, description=""))
    

SeqIO.write(seqs, res_path + "/reference_sequences/aso_targets.fasta", "fasta")
output_df.to_csv(res_path + "/outputs/result_table.tsv", sep="\t", index=False)
#subprocess.run(["Rscript",  "./scripts/melting.R", res_path + "/outputs/result_table.tsv"])




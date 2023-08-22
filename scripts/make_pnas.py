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

# import used gene length:
length = int(sys.argv[1])
res_path = sys.argv[2]
try:
    bases_before = int(sys.argv[3])
except IndexError:
    bases_before = length-3

print(bases_before)

# import fasta file:
target_regions = list(SeqIO.parse(res_path + "/reference_sequences" + "/targetgene_startreg.fasta", "fasta"))

list_aso_sequences = []
list_aso_targets = []
count = 1
sd = "GAGG"

# go through target regions and
for r in range(len(target_regions)):
    gene = re.sub("::.*", "", target_regions[r].id)
    seq = target_regions[r].seq
    seqlist = [i for i in seq.transcribe().__str__()]  # add sequence list for heatmap
    seqlist_pnas = [i for i in seq.complement().__str__()]
    heatmap_list = []
    heatmap_annot = []
    output_df = pd.DataFrame(columns=["ASO", "ASO_seq", "target_seq", "location",
                                      "SC_bases", "pur_perc", "long_pur_stretch", "Tm",
                                      "OT_tot", "OT_TIR"])
    # identify location of SD regions:

    if seq[15:28].find(sd) != -1:
        sd_pos = seq[15:28].find(sd) + 15

    for i in range(30-bases_before, 31):   # 30 is start of cds
        s = seq[i:i+length]

        aso = s.reverse_complement()

        maxcomp = SequenceMatcher(None, aso, s).find_longest_match(0, len(aso), 0, len(s)).size
        
        # design PNA for regions overlapping SC or SD region:
        if maxcomp < length*0.6:
            aso_name = gene+"_ASO_"+str(count).zfill(3)
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
            pur_perc = "{:.2f}".format((pur/len(aso.__str__()))*100)

            # add to dataframe:
            added_row = pd.Series([aso_name, aso.__str__(), s.transcribe().__str__(), str(i-30) + ";" + str(i-30+length),
                                   maxcomp, pur_perc, longest_purine_stretch, None, None, None], index=output_df.columns)
            #output_df = output_df.append(added_row, ignore_index=True)
            output_df = pd.concat([output_df, pd.DataFrame([added_row])], ignore_index=True)
            # save ASO sequences:
            list_aso_sequences += [SeqRecord(aso, id=aso_name, description="")]
            # for fasta with reverse as well to get mismatches with seqmap:
            list_aso_targets += [SeqRecord(s, id=aso_name, description="")]

            heatmap_list += [i*[0] + length*[1] + (len(seq)-i-length)*[0]]
            heatmap_list[count-1][30:33] = [0.4 + i for i in heatmap_list[count-1][30:33]]
            if seq[15:28].find(sd) != -1:
                heatmap_list[count - 1][sd_pos:sd_pos+len(sd)] = [0.2 + i for i in heatmap_list[count - 1][sd_pos:sd_pos+len(sd)] ]
            heatmap_annot += [i*[""] + seqlist_pnas[i:i+length] + (len(seq)-i-length)*[""]]
            count += 1
    heatmap_array = np.array(heatmap_list)

    SeqIO.write(list_aso_sequences, res_path + "/reference_sequences/aso_sequences.fasta", "fasta")
    SeqIO.write(list_aso_targets, res_path + "/reference_sequences/aso_targets.fasta", "fasta")
    # safe output df:
    output_df.to_csv(res_path + "/outputs/result_table.tsv", sep="\t", index=False)
    # run rscript:
    print(res_path + "/outputs/result_table.tsv")
    subprocess.run(["Rscript", "./scripts/melting.R", res_path + "/outputs/result_table.tsv"])

    # seaborn:
    fig, ax = plt.subplots(figsize=(15, count/3))
    ax = sns.heatmap(heatmap_array, linewidths=.2, linecolor="black",
                     annot=np.array(heatmap_annot), fmt="", cmap="Blues",
                     cbar=False)
    ax.set_xticks(np.arange(0.5, len(seqlist), 1))
    ax.set_xticklabels(seqlist, rotation=0, fontsize=12)
    ax.set_xlabel("mRNA (5' to 3')", fontsize=12)
    ylabels = ["ASO_" + str(i + 1).zfill(3) for i in range(count-1)]
    ax.set_yticklabels(ylabels, rotation=0)
    ax.set_ylabel("ASO sequence", fontsize=12)
    plt.tight_layout()
    plt.savefig(res_path + '/outputs/heatmap.png', dpi=500)
    plt.savefig(res_path + '/outputs/heatmap.svg')

    count = 0





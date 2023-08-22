import os
import matplotlib.pyplot as plt
import re
import pandas as pd
import seaborn as sns
import numpy as np
import sys
import six


ot_table = sys.argv[1] + "/offtargets_fulltranscripts_sorted.tab"
all_off_targets = pd.read_table(ot_table, sep='\t', index_col=False)

all_off_targets["ASO"] = all_off_targets["probe_id"].replace("_rev", "", regex=True)

# only use the mismatches with >+7 cons. mm:
all_off_targets = all_off_targets[all_off_targets["longest_stretch"] >= 7]

# add output df used for other things, e.g. Tm:
output_df = pd.read_csv(sys.argv[1] + "/result_table.tsv", sep="\t", index_col=None)

# add variable for whether it is in the TIR (-20 to +5 start)
all_off_targets.loc[all_off_targets["trans_coord"].isin(range(-20, 6)), "TIR"] = "TIR"
all_off_targets.loc[~all_off_targets["trans_coord"].isin(range(-20, 6)), "TIR"] = "not in TIR"

all_off_targets.to_csv(sys.argv[1] + "/offtargets_fulltranscripts_sorted.csv")
all_off_targets[all_off_targets["TIR"] == "TIR"].to_csv(sys.argv[1] + "/offtargets_startregions_sorted.csv")

df_plot = pd.DataFrame(columns=["ASO", "off-target type", "transcripts", "counts", "target sequence"])

# create dataframe:
for i in all_off_targets["ASO"].unique():
    target_seq = all_off_targets[all_off_targets["ASO"] == i].iloc[0, ]["probe_seq"]
    aso_n = re.sub(".*_(ASO.*)", "\\1", i)
    ot_aso = all_off_targets[all_off_targets["ASO"] == i]
    num_tot_ot = ot_aso.shape[0]-1
    num_tir_ot = ot_aso[ot_aso["TIR"] == "TIR"].shape[0]-1
    df_plot = pd.concat([df_plot, pd.DataFrame([[aso_n, "OT in transcriptome", "whole transcriptome", num_tot_ot,
                                                 target_seq]], columns=df_plot.columns)])
    df_plot = pd.concat([df_plot, pd.DataFrame([[aso_n, "OT in TIR regions", "start regions", num_tir_ot,
                                                 target_seq]], columns=df_plot.columns)])
    # change output:
    output_df.loc[ i, "OT_tot"] = num_tot_ot
    output_df.loc[i, "OT_TIR"] = num_tir_ot


# visualize:
def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=(25, data.shape[1]))
        ax.axis('off')

    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in six.iteritems(mpl_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0] % len(row_colors)])
    return ax


def show_on_single_plot(ax, maxlen):
    for p in ax.patches:
        _x = p.get_x() + p.get_width() / 2
        _y = p.get_y() + p.get_height() + maxlen/120
        value = '{:.0f}'.format(p.get_height())
        ax.text(_x, _y, value, ha="center", fontsize=15)


def create_ot_barplot(dataframe, title, filepath):
    # get the number of asos:
    sns.set_style("whitegrid")
    plt.figure(figsize=(16, 11))
    plt.xlabel('xlabel', fontsize=20, fontweight='bold')
    plt.ylabel('ylabel', fontsize=20, fontweight='bold')
    plt.title('titlelabel', fontsize=30, fontweight='bold')

    bp = sns.barplot(x="ASO", y="counts", hue="off-target type", data=dataframe,
                     palette='Blues_r')
    bp.set(xlabel='ASO sequence', ylabel='Number of off-targets')
    show_on_single_plot(bp, maxlen=dataframe["counts"].max())
    bp.set_title(title, fontsize=30, fontweight='bold')
    bp.set_xticklabels(bp.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize=15)
    bp.legend(title="off-target type", fontsize=15, title_fontsize=20)
    bp.figure.savefig(filepath + ".png")
    bp.figure.savefig(filepath + ".svg")
    plt.clf()


create_ot_barplot(df_plot, "Critical off-targets of ASOs", sys.argv[1] + "/plot_ots_whole_transcriptome")


output_df["OT_tot"] = output_df["OT_tot"]+2#.astype(int)
output_df["OT_TIR"] = output_df["OT_TIR"]+2#.astype(int)
ax = render_mpl_table(output_df, header_columns=0, col_width=4.0)
ax.figure.savefig(sys.argv[1] + "/result_table.png", bbox_inches='tight')
ax.figure.savefig(sys.argv[1] + "/result_table.svg", bbox_inches='tight')

# save df for download/visualization:
output_df.to_csv(sys.argv[1] + "/result_table.csv", sep=",")
os.remove(sys.argv[1] + "/result_table.tsv")


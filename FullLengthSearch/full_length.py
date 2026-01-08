from urllib.parse import urlparse
from Bio import SeqIO
from Bio import AlignIO
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import itertools
import os, gzip, shutil, glob, json, ast, sys, tempfile, re
from io import StringIO
from Bio.Blast import NCBIXML

import tools

# 1. foldseek search server 
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# json_folder = os.path.join(BASE_DIR, "foldseekResults")
# json_files = glob.glob(os.path.join(json_folder, "*.json")) 
# df = tools.parse_foldseek_json(json_files)

# 2. remove redundant sequences
# with open("fulllength2199.fasta", "w") as f:
#     for _, row in df[["target", "tSeq"]].iterrows():
#         f.write(f">{row['target']}\n")
#         f.write(f"{row['tSeq']}\n")

# new_df = tools.parse_CDHIT_csv(df=df, fp=os.path.join(BASE_DIR, "fulllength2199_CDHIT.csv"))
# new_df.to_csv("fulllength2199_CDHIT949.csv")
new_df_tmp = pd.read_csv(os.path.join(BASE_DIR, "fulllength2199_CDHIT949.csv"))
print(len(new_df_tmp)) # 949
new_df = new_df_tmp[(new_df_tmp["prob"] >= 0.9) & (new_df_tmp["e_val"] <= 0.01)]
print(len(new_df)) # 879


# 3. distribution of structure alignment hits
mapping = {i: 0 for i in range(2521)}
for _, row in new_df.iterrows():
    qStartPos = row["qStartPos"]-1
    qEndPos = row["qEndPos"]-1
    for i in range(qStartPos, qEndPos+1):
        mapping[i] += 1

# tools.VIZ(mapping=mapping, windows=10, yname="Number of structure-based <br> aligned sequences", xname="Q92508 residue position", svgname=os.path.join(BASE_DIR, "seq_fulllength.svg"), save=True)

# 4. Prob and e value
# rest_df = new_df_tmp[~((new_df_tmp["prob"] >= 0.9) & (new_df_tmp["e_val"] <= 0.01))]
# print(len(rest_df))
# Prob_mapping = {i: [] for i in range(2521)}
# E_mapping = {i: [] for i in range(2521)}
# for _, row in rest_df.iterrows():
#     qStartPos = row["qStartPos"]-1
#     qEndPos = row["qEndPos"]-1
#     prob = row["prob"]
#     e = row["e_val"]
#     for i in range(qStartPos, qEndPos+1):
#         Prob_mapping[i].append(float(prob))
#         E_mapping[i].append(float(e))



# tools.VIZ_prob_e(mapping=Prob_mapping, yname="Prob")
# tools.VIZ_prob_e(mapping=E_mapping, yname="E-value")

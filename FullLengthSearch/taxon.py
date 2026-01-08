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

# with open("foldseek_piezo949.fasta", "w") as f:
#     for _, row in new_df.iterrows():
#         ID = row["target"]
#         seq = df[df["target"]==ID]["tSeq"].tolist()[0]
#         f.write(f">{ID}\n{seq}\n")

# with open("foldseek_piezo949_AF_A0.txt", "w") as f:
#     for _, row in new_df.iterrows():
#         ID = row["target"]
#         if ID[:2] == "AF":
#             name = ID.split("-")[1]
#             f.write(f"{name}\n")
#         elif ID[:2] == "A0":
#             f.write(f"{ID}\n")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
new_df = pd.read_csv(os.path.join(BASE_DIR, "fulllength2199_CDHIT949.csv"))

id_taxid = {}

for record in SeqIO.parse(os.path.join(BASE_DIR, "idmapping_2025_12_26.fasta/idmapping_2025_12_26.fasta"), "fasta"):
    header = record.description
    uid = record.id.split("|")[1]

    for field in header.split():
        if field.startswith("OX="):
            id_taxid[uid] = field.replace("OX=", "")

taxid_list = []
ID_list = []
others_list = []
for _, row in new_df.iterrows():
    ID = row["target"]
    if ID[:2] == "AF":
        name = ID.split("-")[1]
        taxid = id_taxid.get(name, None)
        if taxid is None:
            others_list.append(ID)
        else:
            taxid_list.append(taxid)
            ID_list.append(name)
    elif ID[:2] == "A0":
        taxid = id_taxid.get(ID, None)
        if taxid is None:
            others_list.append(ID)
        else:
            taxid_list.append(taxid)
            ID_list.append(ID)
    else:
        others_list.append(ID)


print(len(others_list), len(ID_list), len(taxid_list))

# with open("foldseek_piezo949_exclude_111.txt", "w") as f:
#     for ID in others_list:
#         seq = df[df["target"]==ID]["tSeq"].tolist()[0]
#         f.write(f">{ID}\n{seq}\n")

data = tools.parse_blastxml2(os.path.join(BASE_DIR, "MYYZBU9W016-Alignment.xml"))

query_tax_dict={}
for item in data:
    query_title = item["query_title"]
    hits = item["hits"]
    if len(hits) == 0:
        # print(query_title, hits)
        # MGYP000926044221 [] 2600159 https://www.ebi.ac.uk/Tools/hmmer/
        # MGYP003976619017 [] 30611 https://www.ebi.ac.uk/Tools/hmmer/
        # MGYP000968413431 [] 9694 https://www.ebi.ac.uk/Tools/hmmer/
        # A0A345MQM9 [] 322159 in blastp
        # MGYP003539167077 [] 317549 https://www.ebi.ac.uk/Tools/hmmer/
        if query_title == "MGYP000926044221":
            query_tax_dict[query_title] = 2600159
        if query_title == "MGYP003976619017":
            query_tax_dict[query_title] = 30611
        if query_title == "MGYP000968413431":
            query_tax_dict[query_title] = 9694
        if query_title == "A0A345MQM9":
            query_tax_dict[query_title] = 322159
        if query_title == "MGYP003539167077":
            query_tax_dict[query_title] = 317549
    else:
        # print(query_title, hits[0]["taxid"])
        query_tax_dict[query_title] = hits[0]["taxid"]

    # query_tax_dict["AF-Q9H5I5-2-F1-model_v6"] = 9606 # uniprot website
    # query_tax_dict["AF-A0AB32U1Q1-F1-model_v6"] = 7955 # uniprot website
    # query_tax_dict["AF-Q4DXU9-F1-model_v6"] = 353153 # uniprot website
for i, j in zip(ID_list, taxid_list):
    query_tax_dict[i] = j





# get all taxID
# print(query_tax_dict)
print(len(query_tax_dict))
new_df_tmp = pd.read_csv(os.path.join(BASE_DIR, "fulllength2199_CDHIT949.csv"))
print(len(new_df_tmp)) # 949
tax_df_tmp = new_df_tmp[(new_df_tmp["prob"] >= 0.9) & (new_df_tmp["e_val"] <= 0.01)]
print(len(tax_df_tmp)) # 879
# print(tax_df_tmp)
taxID_list = []
query_list = []
for _, row in tax_df_tmp.iterrows():
    q = row["target"]
    try:
        taxID = query_tax_dict[q]
    except KeyError:
        q_split = q.split("-")[1]
        taxID = query_tax_dict[q_split]
    taxID_list.append(taxID)
    query_list.append(q)


tax_df = pd.DataFrame({
    "target": query_list,
    "taxID": taxID_list
})

tax_df.to_csv(os.path.join(BASE_DIR, "tax_df_879.csv"))
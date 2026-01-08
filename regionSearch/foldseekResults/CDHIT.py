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
import xml.etree.ElementTree as ET

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(BASE_DIR, "FoldseekSearch_4regions_result.csv"))

df_blade = df[df["queries_header"]=="job_A blade"]
df_beam = df[df["queries_header"]=="job_A beam"]
df_cap = df[df["queries_header"]=="job_A cap"]
df_pore = df[df["queries_header"]=="job_A pore"]

with open("FoldseekSearch_4regions_result_blade.fasta", "w") as f:
    for _, row in df_blade[["target", "tSeq"]].iterrows():
        f.write(f">{row['target']}\n")
        f.write(f"{row['tSeq']}\n")

with open("FoldseekSearch_4regions_result_beam.fasta", "w") as f:
    for _, row in df_beam[["target", "tSeq"]].iterrows():
        f.write(f">{row['target']}\n")
        f.write(f"{row['tSeq']}\n")
        
with open("FoldseekSearch_4regions_result_cap.fasta", "w") as f:
    for _, row in df_cap[["target", "tSeq"]].iterrows():
        f.write(f">{row['target']}\n")
        f.write(f"{row['tSeq']}\n")

with open("FoldseekSearch_4regions_result_pore.fasta", "w") as f:
    for _, row in df_pore[["target", "tSeq"]].iterrows():
        f.write(f">{row['target']}\n")
        f.write(f"{row['tSeq']}\n")
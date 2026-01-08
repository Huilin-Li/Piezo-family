from pycdhit import read_fasta, cd_hit, read_clstr, CDHIT


df_in = read_fasta("./up2now_piezo55.fasta")
cdhit = CDHIT(prog="cd-hit", path="~/tools/cd-hit-v4.8.1-2019-0228")
df_out, df_clstr = cdhit.set_options(c=0.9).cluster(df_in)
print(df_clstr)

df_clstr.to_csv("./up2now_piezo_distinct_CDHIT.csv")
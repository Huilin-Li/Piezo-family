import pandas as pd
from urllib.parse import urlparse
from Bio import SeqIO
from Bio import AlignIO
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import json
from plotly.subplots import make_subplots

def parse_3line(fp):
    with open(fp) as f:
        lines = f.readlines()
    seq = lines[1][:-1]
    IMO = lines[2][:-1]
    df = pd.DataFrame({
        "pos": list(range(len(seq))),
        "aa":  list(seq),
        "IMO": list(IMO)
    })
    return df

def viz_TM(tmdf, xaxis_title):
    x_pos=0
    n = 60
    figTMHMMeach = go.Figure()
    
    for x in range(0, len(tmdf)+1, 50):
        figTMHMMeach.add_trace(go.Scatter(
            x=[x, x],
            y=[-1.1, 1.1],
            mode="lines",
            line=dict(color="rgba(128,128,128,0.5)", width=1),
            showlegend=False,
            hoverinfo="skip"
        ))
    
    # ---- underlay horizontal lines as traces ----
    for y in [1.1, 1, 0, -1, -1.1]:
        figTMHMMeach.add_trace(go.Scatter(
            x=[0, len(tmdf)-1],
            y=[y, y],
            mode="lines",
            line=dict(color="rgba(128,128,128,0.5)", width=1),
            showlegend=False,
            hoverinfo="skip"
        ))
    
    x_positions = list(range(len(tmdf)))

    for x_pos, imo in zip(x_positions, tmdf["IMO"]):
        if imo == "M":
            figTMHMMeach.add_trace(go.Scatter(
                x=[x_pos] * n,
                y=np.linspace(-1, 1, n),
                mode="markers",
                marker=dict(color="red"),
                showlegend=False
            ))

    
    for x_pos, imo in zip(x_positions, tmdf["IMO"]):
        if imo == "I":
            figTMHMMeach.add_trace(go.Scatter(
                x=[x_pos],
                y=[-1],
                mode="markers",
                marker=dict(color="pink", size=7),
                showlegend=False
            ))

    for x_pos, imo in zip(x_positions, tmdf["IMO"]):
        if imo == "O":
            figTMHMMeach.add_trace(go.Scatter(
                x=[x_pos],
                y=[1],
                mode="markers",
                marker=dict(color="blue", size=7),
                showlegend=False
            ))
    
    figTMHMMeach.update_xaxes(range=[-0.2, len(tmdf)+1], title_text=xaxis_title)
    figTMHMMeach.update_yaxes(
        range=[-1.1, 1.1],
        tickmode="array",
        tickvals=[-1, 0, 1],
        ticktext=["I", "M", "O"],
        showticklabels=True,
        ticks="outside",
        showline=False,
    )
    figTMHMMeach.update_layout(
        height=200,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=0, r=0, t=20, b=20)
    )
    
    # figTMHMMeach.update_yaxes(
    #     range=[-1.1, 1.1],
    #     anchor="free",
    #     position=0,
    #     side="left",
    #     tickmode="array",
    #     tickvals=[-1, 0, 1],
    #     ticktext=["I", "M", "O"],
    #     showticklabels=True,
    #     ticks="outside",
    #     showline=False,
    # )
    # figTMHMMeach.update_layout(
    #     xaxis_title=xaxis_title,
    #     yaxis_title="",
    #     xaxis_title_font=dict(size=11),
    #     yaxis_title_font=dict(size=11),
    #     width=500,
    #     height=200,
    #     paper_bgcolor="rgba(0,0,0,0)",
    #     plot_bgcolor="rgba(0,0,0,0)",
    #     margin=dict(l=20, r=20, t=20, b=20)
    #     )
    # figTMHMMeach.write_image(name)
    # figTMHMMeach.show()
    return figTMHMMeach


A0A4_tm_fp = "E:/CLUSTERwork/PIEZO1_EVO/blade/DeepTMHMM/A0A48deeptmhmm/predicted_topologies.3line"
mgyp_tm_fp = "E:/CLUSTERwork/PIEZO1_EVO/blade/DeepTMHMM/MGYPdeeptmhmm/predicted_topologies.3line"
blade_tm_fp = "E:/CLUSTERwork/PIEZO1_EVO/blade/DeepTMHMM/blade_deeptmhmm/predicted_topologies.3line"

A0A4_tm_df = parse_3line(A0A4_tm_fp)
mgyp_tm_df = parse_3line(mgyp_tm_fp)
blade_tm_df = parse_3line(blade_tm_fp)



fig1 = viz_TM(blade_tm_df, "hPiezo1 blade region residue position")
fig2 = viz_TM(mgyp_tm_df, "MGYP000926044221 residue position")
fig3 = viz_TM(A0A4_tm_df, "A0A482RRZ8 residue position")

big = make_subplots(
    rows=3, cols=1,
    shared_xaxes=False,   
    vertical_spacing=0.15,

)


for tr in fig1.data: big.add_trace(tr, row=1, col=1)
for tr in fig2.data: big.add_trace(tr, row=2, col=1)
for tr in fig3.data: big.add_trace(tr, row=3, col=1)


big.update_xaxes(range=fig1.layout.xaxis.range, title_text=fig1.layout.xaxis.title.text, row=1, col=1)
big.update_yaxes(range=fig1.layout.yaxis.range, tickmode="array", tickvals=[-1,0,1], ticktext=["I","M","O"], row=1, col=1)

big.update_xaxes(range=fig2.layout.xaxis.range, title_text=fig2.layout.xaxis.title.text, row=2, col=1)
big.update_yaxes(range=fig2.layout.yaxis.range, tickmode="array", tickvals=[-1,0,1], ticktext=["I","M","O"], row=2, col=1)

big.update_xaxes(range=fig3.layout.xaxis.range, title_text=fig3.layout.xaxis.title.text, row=3, col=1)
big.update_yaxes(range=fig3.layout.yaxis.range, tickmode="array", tickvals=[-1,0,1], ticktext=["I","M","O"], row=3, col=1)

big.update_layout(
    height=3*220,
    width=500,  
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    margin=dict(l=0, r=0, t=20, b=10)
)
big.write_image("E:/CLUSTERwork/PIEZO1_EVO/blade/DeepTMHMM/TMVIZ.svg")
big.show()
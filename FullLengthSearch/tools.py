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

def convert_hrefs(href):
    df = pd.DataFrame(href)
    df["db_index"] = df.groupby("label").cumcount() + 1

    row_data = {}
    for _, row in df.iterrows():
        label = row["label"]
        idx = row["db_index"]
        row_data[f"{label}{idx}"] = 1              
        row_data[f"{label}link{idx}"] = row["href"] 

    return row_data


def parse_foldseek_json(json_files):
    ROWSforDF = []
    for json_fp in json_files:
        with open(json_fp) as f:
            json_d = json.load(f)
            #print(json_d[0].keys())
            queries_header = json_d[0]["queries"][0]["header"]
            queries_seq = json_d[0]["queries"][0]["sequence"]
            mode = json_d[0]["mode"]
            # print(json_fp,"###",queries_header, "###",mode)
            print("###",queries_header, "###",mode)
            res = json_d[0]["results"]
            for i in range(9):
                db = res[i]["db"]
                alignments_count = len(res[i]["alignments"])
                # print(db, len(res[i]["alignments"]))
                if alignments_count == 0:
                    pass
                else:
                    for j in range(alignments_count):
                        target = res[i]["alignments"][str(j)][0]["target"]
                        targte_dict = {"target": target}
                        prob = res[i]["alignments"][str(j)][0]["prob"]
                        prob_dict = {"prob": prob}
                        e_val = res[i]["alignments"][str(j)][0]["eval"]
                        e_val_dict = {"e_val": e_val}
                        score = res[i]["alignments"][str(j)][0]["score"]
                        score_dict = {"score": score}
                        qStartPos = res[i]["alignments"][str(j)][0]["qStartPos"]
                        qStartPos_dict = {"qStartPos": qStartPos}
                        qEndPos = res[i]["alignments"][str(j)][0]["qEndPos"]
                        qEndPos_dict = {"qEndPos": qEndPos}
                        qLen = res[i]["alignments"][str(j)][0]["qLen"]
                        qLen_dict = {"qLen": qLen}
                        dbLen = res[i]["alignments"][str(j)][0]["dbLen"]
                        dbLen_dict = {"dbLen": dbLen}
                        qAln = res[i]["alignments"][str(j)][0]["qAln"]
                        qAln_dict = {"qAln": qAln}
                        dbAln = res[i]["alignments"][str(j)][0]["dbAln"]
                        dbAln_dict = {"dbAln": dbAln}
                        tSeq = res[i]["alignments"][str(j)][0]["tSeq"]
                        tSeq_dict = {"tSeq": tSeq}
                        taxId = res[i]["alignments"][str(j)][0].get("taxId", "None")
                        taxId_dict = {"taxId": taxId}
                        taxName = res[i]["alignments"][str(j)][0].get("taxName", "None")
                        taxName_dict = {"taxName": taxName}
                        description = res[i]["alignments"][str(j)][0].get("description", "None")
                        description_dict = {"description": description}
                        href0 = res[i]["alignments"][str(j)][0]["href"]
                        
                        queries_header_dict = {"queries_header": queries_header}
                        queries_seq_dict = {"queries_seq": queries_seq}
                        mode_dict = {"mode": mode}
                        db_dict = {"db": db}

                        if type(href0) != list:
                            # {'gmgc.embl.de', 'esmatlas.com', 'modelarchive.org', 'www.rcsb.org', 'www.cathdb.info', 'predictomes.org'}
                            u = urlparse(href0)
                            netloc = u.netloc
                            href_dict = {netloc: href0}
                        else:
                            href_dict = convert_hrefs(href0)
                            
                        dict_list = [queries_header_dict, queries_seq_dict, mode_dict, db_dict, targte_dict, prob_dict, e_val_dict, score_dict, qStartPos_dict, qEndPos_dict, qLen_dict, dbLen_dict, qAln_dict, dbAln_dict, tSeq_dict, taxId_dict, taxName_dict, description_dict, href_dict]

                        oneRow = {}
                        for d in dict_list:
                            oneRow |= d 

                        ROWSforDF.append(oneRow)
    df = pd.DataFrame(ROWSforDF)
    return df

def parse_CDHIT_csv(df, fp):
    rows = []
    after_CDHIT_df = pd.read_csv(fp)
    for cluster in set(after_CDHIT_df["cluster"].tolist()):
        cluster_df = after_CDHIT_df[after_CDHIT_df["cluster"]==cluster].sort_values(by="size", ascending=False)
        identifier = cluster_df["identifier"].tolist()[0]
        filtered_df = df[df["target"].str.contains(identifier)].head(1)

        for _, row in filtered_df.iterrows():
            rows.append({
                "target": row["target"],
                "prob": row["prob"],
                "e_val": row["e_val"],
                "qStartPos": row["qStartPos"],
                "qEndPos": row["qEndPos"]
            })

    new_df = pd.DataFrame(rows)
    return new_df



def VIZ(mapping, windows, yname, xname, svgname, save=True):
    BladeTHUs1_7 = [1,1304]
    Beam = [1305,1367]
    Central_plug = [1368,1401]
    Lateral_plug=[1402, 1411]
    Latch=[1412,1425]
    Latch2Clasp=[1426,1514]
    Clasp=[1515,1570]
    Clasp2THU8=[1571,1656]
    BladeTHUs8_9=[1657,2098]
    Anchor=[2099,2165]
    OuterHelix=[2166,2193]
    Cap=[2194,2439]
    InnerHelix=[2440,2473]
    CTD=[2474,2521]
    
    colors_12 = px.colors.qualitative.Light24[:14]
    
    regions = [
        ("BladeTHUs1_7",   BladeTHUs1_7, colors_12[0]),
        ("Beam",           Beam, colors_12[1]),
        ("Central_plug",   Central_plug, colors_12[2]),
        ("Lateral_plug",   Lateral_plug, colors_12[3]),
        ("Latch",          Latch, colors_12[4]),
        ("Latch2Clasp",    Latch2Clasp, colors_12[5]),
        ("Clasp",          Clasp, colors_12[6]),
        ("Clasp2THU8",     Clasp2THU8, colors_12[7]),
        ("BladeTHUs8_9",   BladeTHUs8_9, colors_12[8]),
        ("Anchor",         Anchor, colors_12[9]),
        ("OuterHelix",     OuterHelix, colors_12[10]),
        ("Cap",            Cap, colors_12[11]),
        ("InnerHelix",     InnerHelix, colors_12[12]),
        ("CTD",            CTD, colors_12[13]),
    ]
    x = np.arange(1, 2522)
    colors = np.full(len(x), "lightgray", dtype=object)
    
    for _, start_end, color in regions:
        mask = (x >= start_end[0]) & (x <= start_end[-1])
        colors[mask] = color
        
    y = list(mapping.values())
    y_arr = np.array(y)
    
    # print("len(colors), len(y_arr), len(x)", len(colors), len(y_arr), len(x))

    vz_df = pd.DataFrame({
        'x': x,
        'y_arr': y_arr,
        'color': colors
        })

    # regression
    from sklearn.tree import DecisionTreeRegressor
    from sklearn.ensemble import AdaBoostRegressor
    regr = AdaBoostRegressor(
        DecisionTreeRegressor(max_depth=4), n_estimators=300, random_state=42
    )
    regr.fit(x.reshape(-1, 1), y_arr)
    y_fit = regr.predict(x.reshape(-1, 1))
    fit_df = pd.DataFrame({
        'x': x,
        'y_fit': y_fit,
        'color': colors
    })
    window_size = windows  # Adjust the window size as needed
    y_smooth = np.convolve(y_arr, np.ones(window_size)/window_size, mode='valid')

    figBar = px.bar(vz_df, y='y_arr', x='x',  color="color", color_discrete_map="identity")
    figBar.update_layout(
        xaxis_title=xname,
        yaxis_title=yname,
        bargap=0,
        bargroupgap=0,
        width=1200,
        height=300,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=20, r=20, t=0, b=0),
        showlegend=True
    )
    figBar.add_trace(
        go.Scatter(
            x=vz_df['x'][window_size-1:],  # Align x-axis to the smoothed data (to match length)
            y=y_smooth,
            mode='lines',
            name='Smoothed Trend',
            line=dict(color='black', width=0.7),
            showlegend=False
        )
    )
    for region_name, start_end, color in regions:
        figBar.add_trace(
            go.Scatter(
                x=start_end,  # Position the region on the x-axis
                y=[max(vz_df['y_arr'])] * 2,  # Use max y to position the legend
                mode='lines',
                line=dict(color=color, width=10),
                name=region_name
            )
        )
    figBar.add_vline(x=1305, line_width=2, line_dash="dash", line_color="black")
    figBar.add_vline(x=2156, line_width=2, line_dash="dash", line_color="black")
    figBar.add_vline(x=2210, line_width=2, line_dash="dash", line_color="black")
    figBar.add_vline(x=2428, line_width=2, line_dash="dash", line_color="black")

    # Update layout to add the legend at the bottom
    figBar.update_layout(
        xaxis_title=xname,
        yaxis_title=yname,
        xaxis_title_font=dict(size=11),
        yaxis_title_font=dict(size=11),
        bargap=0,
        bargroupgap=0,
        width=1100,
        height=300,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=20, r=20, t=0, b=100),  # Increased bottom margin to fit the legend
        showlegend=True,
        legend=dict(
            orientation='h',  # Horizontal legend at the bottom
            x=0,              # Align the legend to the left
            y=-0.5,           # Place the legend below the plot
            xanchor='left',
            yanchor='top',
            font=dict(size=10)
        )
    )
    figBar.update_xaxes(
        tickmode="array",
        tickvals=[500, 1000, 1304, 1500, 2000, 2156, 2210, 2428, 2500],              # only one tick
        ticktext=["500", "1000", "<b>1305</b>", "1500", "2000", "<b>2156</b>", "<b>2210</b>", "<b>2428</b>", "2500"],     # bold number (HTML works)
        tickfont=dict(
            color="black",
            size=10
        ),
        showline=True,
        linecolor="black",
        linewidth=1,
        ticks="outside"
    )

    if save:
        figBar.write_image(svgname)
        #figBar.write_html(svgname)
    figBar.show()

def parse_blastxml2(path: str):
    NS = {"n": "http://www.ncbi.nlm.nih.gov"} 
    tree = ET.parse(path)
    root = tree.getroot()

    records = []

    # Find each Search block
    for search in root.findall(".//n:Results/n:search/n:Search", NS):
        rec = {
            # "query_id": search.findtext("n:query-id", default="", namespaces=NS),
            "query_title": search.findtext("n:query-title", default="", namespaces=NS),
            #"query_len": int(search.findtext("n:query-len", default="0", namespaces=NS)),
            "hits": []
        }

        for hit in search.findall("n:hits/n:Hit", NS):
            descr = hit.find("n:description/n:HitDescr", NS)

            taxid_text = descr.findtext("n:taxid", default="", namespaces=NS) if descr is not None else ""
            taxid = int(taxid_text) if taxid_text.isdigit() else None

            hit_obj = {
                #"num": int(hit.findtext("n:num", default="0", namespaces=NS)),
                #"accession": descr.findtext("n:accession", default="", namespaces=NS) if descr is not None else "",
                #"id": descr.findtext("n:id", default="", namespaces=NS) if descr is not None else "",
                #"title": descr.findtext("n:title", default="", namespaces=NS) if descr is not None else "",
                "taxid": taxid,
                "sciname": descr.findtext("n:sciname", default="", namespaces=NS) if descr is not None else "",
                #"len": int(hit.findtext("n:len", default="0", namespaces=NS)),
                "hsps": []
            }

            for hsp in hit.findall("n:hsps/n:Hsp", NS):
                hit_obj["hsps"].append({
                    #"num": int(hsp.findtext("n:num", default="0", namespaces=NS)),
                    "evalue": float(hsp.findtext("n:evalue", default="0", namespaces=NS)),
                    #"bit_score": float(hsp.findtext("n:bit-score", default="0", namespaces=NS)),
                    # "score": int(hsp.findtext("n:score", default="0", namespaces=NS)),
                    # "identity": int(hsp.findtext("n:identity", default="0", namespaces=NS)),
                    # "align_len": int(hsp.findtext("n:align-len", default="0", namespaces=NS)),
                    # "query_from": int(hsp.findtext("n:query-from", default="0", namespaces=NS)),
                    # "query_to": int(hsp.findtext("n:query-to", default="0", namespaces=NS)),
                    # "hit_from": int(hsp.findtext("n:hit-from", default="0", namespaces=NS)),
                    # "hit_to": int(hsp.findtext("n:hit-to", default="0", namespaces=NS)),
                    # optional huge strings:
                    # "qseq": hsp.findtext("n:qseq", default="", namespaces=NS),
                    # "hseq": hsp.findtext("n:hseq", default="", namespaces=NS),
                })

            rec["hits"].append(hit_obj)

        records.append(rec)

    return records

def VIZ_prob_e(mapping, yname, scatter=True, avgline=False):
    BladeTHUs1_7 = [1,1304]
    Beam = [1305,1367]
    Central_plug = [1368,1401]
    Lateral_plug=[1402, 1411]
    Latch=[1412,1425]
    Latch2Clasp=[1426,1514]
    Clasp=[1515,1570]
    Clasp2THU8=[1571,1656]
    BladeTHUs8_9=[1657,2098]
    Anchor=[2099,2165]
    OuterHelix=[2166,2193]
    Cap=[2194,2439]
    InnerHelix=[2440,2473]
    CTD=[2474,2521]
    
    colors_12 = px.colors.qualitative.Light24[:14]
    
    regions = [
        ("BladeTHUs1_7",   BladeTHUs1_7, colors_12[0]),
        ("Beam",           Beam, colors_12[1]),
        ("Central_plug",   Central_plug, colors_12[2]),
        ("Lateral_plug",   Lateral_plug, colors_12[3]),
        ("Latch",          Latch, colors_12[4]),
        ("Latch2Clasp",    Latch2Clasp, colors_12[5]),
        ("Clasp",          Clasp, colors_12[6]),
        ("Clasp2THU8",     Clasp2THU8, colors_12[7]),
        ("BladeTHUs8_9",   BladeTHUs8_9, colors_12[8]),
        ("Anchor",         Anchor, colors_12[9]),
        ("OuterHelix",     OuterHelix, colors_12[10]),
        ("Cap",            Cap, colors_12[11]),
        ("InnerHelix",     InnerHelix, colors_12[12]),
        ("CTD",            CTD, colors_12[13]),
    ]
    x = np.arange(1, 2522)
    colors = np.full(len(x), "lightgray", dtype=object)

    for _, start_end, color in regions:
        mask = (x >= start_end[0]) & (x <= start_end[-1])
        colors[mask] = color

    x_vals = []
    y_vals = []
    color_vals = []

    for i, probs in mapping.items():
        for p in probs:
            x_vals.append(i)
            y_vals.append(p)
            color_vals.append(colors[i])
    vz_df = pd.DataFrame({
        'x': x_vals,
        'y': y_vals,
        'color': color_vals
        })
    if scatter:
        fig = go.Figure(
            go.Scatter(
                x=vz_df["x"],
                y=vz_df["y"],
                mode="markers",
                marker=dict(color=vz_df["color"], size=3, opacity=0.6)
            )
        )
    if avgline:
        line_x_vals = []
        line_y_vals = []
        line_color_vals = []

        for i, v in mapping.items():
            line_x_vals.append(i)
            line_y_vals.append(sum(v) / len(v))
            line_color_vals.append(colors[i])
        line_vz_df = pd.DataFrame({
            'x': line_x_vals,
            'y': line_y_vals,
            'color': line_color_vals
            })
        
        fig = go.Figure(
            go.Scatter(
                x=line_vz_df["x"],
                y=line_vz_df["y"],
                mode="markers+lines",
                marker=dict(color=line_vz_df["color"], size=3, opacity=0.6)
            )
        )


    fig.add_vline(x=1305, line_width=2, line_dash="dash", line_color="black")
    fig.add_vline(x=2156, line_width=2, line_dash="dash", line_color="black")
    fig.add_vline(x=2210, line_width=2, line_dash="dash", line_color="black")
    fig.add_vline(x=2428, line_width=2, line_dash="dash", line_color="black")

    # Update layout to add the legend at the bottom
    fig.update_layout(
        xaxis_title="Q92508 residue position",
        yaxis_title=yname,
        xaxis_title_font=dict(size=11),
        yaxis_title_font=dict(size=11),
        bargap=0,
        bargroupgap=0,
        width=1100,
        height=300,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=20, r=20, t=0, b=100),  # Increased bottom margin to fit the legend
        showlegend=True,
        legend=dict(
            orientation='h',  # Horizontal legend at the bottom
            x=0,              # Align the legend to the left
            y=-0.5,           # Place the legend below the plot
            xanchor='left',
            yanchor='top',
            font=dict(size=10)
        )
    )
    fig.update_xaxes(
        tickmode="array",
        tickvals=[500, 1000, 1304, 1500, 2000, 2156, 2210, 2428, 2500],              # only one tick
        ticktext=["500", "1000", "<b>1305</b>", "1500", "2000", "<b>2156</b>", "<b>2210</b>", "<b>2428</b>", "2500"],     # bold number (HTML works)
        tickfont=dict(
            color="black",
            size=10
        ),
        showline=True,
        linecolor="black",
        linewidth=1,
        ticks="outside"
    )
    fig.show()
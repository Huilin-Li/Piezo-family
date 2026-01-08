import numpy as np
import plotly.graph_objects as go
import ast
from plotly.subplots import make_subplots



def plain_viz(y0, y1, df):
    fig = go.Figure()
    centers_list_tmp = df["center"].tolist()
    # centers_list = [ast.literal_eval(item) for item in centers_list_tmp]
    CENTERS_arr = np.asarray(centers_list_tmp)
    xs = CENTERS_arr[:, 0]
    ys = CENTERS_arr[:, 1]

    fig.add_trace(
            go.Scatter(
                    x=xs,
                    y=ys,
                    mode="lines+markers+text",
                    line=dict(shape="spline", width=2, color="cornflowerblue"),
                    marker=dict(size=10, color="cornflowerblue"),
                    hovertext=df["position"].astype(str) + df["AA"],
                    hoverinfo="text",
                    showlegend=False,
                )
                )
    fig.add_hrect(y0=y0, y1=y1, line_width=0, fillcolor="tan", opacity=0.5, layer="below")
    NtermAA_x, NtermAA_y = xs[0], ys[0]
    fig.add_annotation(x=NtermAA_x, y=NtermAA_y, text="N<sub>2</sub>H-", showarrow=False, xshift=-30)
    CtermAA_x, CtermAA_y = xs[-1], ys[-1]
    fig.add_annotation(x=CtermAA_x, y=CtermAA_y, text="-COOH", showarrow=False, xshift=30)

    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)", 
        plot_bgcolor="rgba(0,0,0,0)", 
        autosize=True,
        margin=dict(l=20, r=20, t=0, b=0),
        showlegend=True
            )
    
    fig.write_html("./topology.html")
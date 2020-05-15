import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.figure_factory as ff
from plotly.offline import iplot
import cufflinks
import pulsar_producer

def read_data(filename):
    df = pd.read_csv((filename+'.txt'), sep='\t', index_col=None)
    return df



vars = pulsar_producer.pulsar_arg_names
del vars [-1]
del vars[11]
del vars[4]
del vars[0]
vars.append("chi")

run = read_data("Results/AllVars100")
run.columns = vars

corrs = run.corr(method = 'spearman')
# figure = ff.create_scatterplotmatrix(
#     run
# )


figure = ff.create_annotated_heatmap(
    z = corrs.values,
    x = list(corrs.columns),
    y = list(corrs.index),
    annotation_text=corrs.round(2).values,
    showscale=True,
    colorscale = 'temps')
figure.show()
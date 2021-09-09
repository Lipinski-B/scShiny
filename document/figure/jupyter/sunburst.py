#!usr/bin/env python3

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import seaborn as sns
import numpy as np

import scipy
import scipy.cluster.hierarchy as hac
from scipy.cluster.hierarchy import dendrogram, linkage

from matplotlib import pyplot as plt

#Make data
def make_metadata(matrix):
    data = pd.read_table(matrix, sep = ',', header = 0)
    data = pd.DataFrame(data, columns=['orig.ident', 'SingleR.calls', 'Phase']).dropna()
    return data
    
#Sunburst
def make_sunburst(data):
    fig = px.sunburst(data, path=data.columns, color='SingleR.calls', 
                      branchvalues='total',template='ggplot2',title='Sunburst Metadata',
                      color_discrete_sequence=px.colors.qualitative.Pastel)

    fig.update_traces(textinfo='label+percent entry')
    fig.update_layout(margin=dict(t=50,l=0,r=0,b=50))
    
    return fig
    
    
data = make_metadata('/home/boris/Documents/analyse/metadata_matrix_FL140304.csv')
make_sunburst(data).savefig('/home/boris/Documents/analyse/sunburst.png')

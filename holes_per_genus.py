import pandas as pd
import numpy as np
from scipy.spatial.distance import hamming
import gudhi as gd
import seaborn as sns
from statistics import mode
import multiprocessing
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import gudhi as gd
# from scipy.spatial.distance import hamming
import plotly.graph_objs as go
import networkx as nx
import plotly.graph_objects as go
import plotly.io as pio
from networkx.utils import not_implemented_for, pairwise
from concurrent.futures import ThreadPoolExecutor, as_completed
import plotly
import plotly.express as px
from module_functions import *



info=pd.read_csv('data/pscplos1246.blast',sep='\t',header=None)
minimum_score=100
info=info[info[11]>=minimum_score]
info.reset_index(drop=True,inplace=True)
info=info.loc[:,0:1]
names=pd.read_csv('data/Actinos.ids',sep='\t',dtype='object',header=None)
info_by_pathway_and_genomes_directions=get_info_by_pathway_and_all_genomes_directions(info)
df_by_pathway,df_by_pathway_drop_duplicate,representative_genomes=get_df_by_pathway(info_by_pathway_and_genomes_directions)
df_by_genus_pathway,df_by_genus_pathway_drop_duplicate,representative_genomes_genus_pathway=get_df_by_genus_pathway(df_by_pathway,names)
complex_genus_pathways=get_complex_by_pathways(df_by_genus_pathway_drop_duplicate)
holes=get_holes_by_pathways(complex_genus_pathways,df_by_genus_pathway_drop_duplicate)


f = open("los1246.txt", "w")
f.write(str(holes))
f.close()



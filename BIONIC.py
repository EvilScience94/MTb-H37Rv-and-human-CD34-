# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 12:44:14 2022

@author: K Praveena
"""

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt

# import BIONIC features
features = pd.read_csv("Mtb_feat64.csv", index_col=0)

# first compute pairwise cosine similarities between genes
features_sim = pd.DataFrame(cosine_similarity(features.values), index=features.index, columns=features.index)

# sparsify the complete graph by removing edges below 0.5 similarity
features_sim[features_sim < 0.5] = 0

# create a `networkx` object
feat_net = nx.from_pandas_adjacency(features_sim)

Bionic = nx.to_pandas_adjacency(feat_net)

Bionic = pd.DataFrame(Bionic)
Bionic.to_csv('Bionic.csv')

Bionic_edge = nx.to_pandas_edgelist(feat_net)
Bionic_edge = pd.DataFrame(Bionic_edge)
Bionic_edge.to_csv('Bionic_edge.csv')

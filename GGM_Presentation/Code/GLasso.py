# %%
import numpy as np 
import matplotlib.pyplot as plt 
from sklearn import covariance 
import pandas as pd 
import networkx as nx
from matplotlib.ticker import LinearLocator
# %%
G = nx.erdos_renyi_graph(8, 0.4, seed=2022)
nx.draw(G)

import numpy as np 
import gudhi as gd
import networkx as nx
import plotly.graph_objects as go
import math
import matplotlib as mpl
import matplotlib
import plotly.io as pio
import matplotlib.pyplot as plt
from GeneralisedFormanRicci.frc import gen_graph
from scipy.sparse import *
from scipy import *

def faces(simplices):
    faceset = set()
    for simplex in simplices:
        numnodes = len(simplex)
        for r in range(numnodes, 0, -1):
            for face in combinations(simplex, r):
                faceset.add(tuple(sorted(face)))
    return faceset

def n_faces(face_set, n):
    return filter(lambda face: len(face)==n+1, face_set)

def boundary_operator(face_set, i):
    source_simplices = list(n_faces(face_set, i))
    target_simplices = list(n_faces(face_set, i-1))
    #print(source_simplices, target_simplices)

    if len(target_simplices)==0:
        S = dok_matrix((1, len(source_simplices)), dtype=np.float64)
        S[0, 0:len(source_simplices)] = 1
    else:
        source_simplices_dict = {source_simplices[j]: j for j in range(len(source_simplices))}
        target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}

        S = dok_matrix((len(target_simplices), len(source_simplices)), dtype=np.float64)
        for source_simplex in source_simplices:
            for a in range(len(source_simplex)):
                target_simplex = source_simplex[:a]+source_simplex[(a+1):]
                i = target_simplices_dict[target_simplex]
                j = source_simplices_dict[source_simplex]
                S[i, j] = -1 if a % 2==1 else 1
    
    return S

data = np.load("guanine.npz", allow_pickle=True)
for i in data["PRO"]:
    coords = i["pos"]

rc = gd.RipsComplex(coords, max_edge_length=1.2)
simplex_tree = rc.create_simplex_tree(max_dimension=2)
val = simplex_tree.get_filtration()
simplices = set()
for v in val:
    simplices.add(tuple(v[0]))


laplacian = np.matmul(boundary_operator(simplices, 2).toarray(), np.transpose(boundary_operator(simplices, 2).toarray()))+np.matmul(np.transpose(boundary_operator(simplices, 1).toarray()), boundary_operator(simplices, 1).toarray())
eigval, eigvec = np.linalg.eigh(laplacian)
print(eigval)

"""
# Below is for visualizing eigenvectors in GEPHI software
G = gen_graph(list(n_faces(simplices, 1)), coords, labels=dict()) # Get the Graph Network of Simplicial Complex
edges = list(n_faces(simplices, 1))
edge_dict = dict()
for i in range(len(eigvec[:, 5])):
    edge_dict[edges[i]] = eigvec[i,5]

H = nx.Graph()
for nn in G.nodes():
    H.add_node(nn, viz={'position': {'x': G.nodes[nn]["coords"][0], 'y':G.nodes[nn]["coords"][1], 'z':G.nodes[nn]["coords"][2]}}) # Following XML format in GEPHI GEXF graphs
for ee in G.edges():
    #print(ee)
    val = abs(edge_dict[ee])
    H.add_edge(ee[0], ee[1], eigvec=abs(edge_dict[ee]), Weight=abs(edge_dict[ee]))

nx.write_gexf(H, "1AXC_CA_eig5_10.0.gexf")
"""


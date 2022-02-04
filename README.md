# Hodge-Theory

This manual is for the code implementation of paper "Hodge Theory for Biomolecular Data Analysis".

# Code Requirements
---
        Platform: Python>=3.6, MATLAB 2016B, GEPHI 0.9.2
        Python Packages needed: math, numpy>=1.19.5, scipy>=1.4.1, scikit-learn>=0.20.3, GUDHI 3.0.0, NetworkX 2.4

# Spectral Eigenvectors 

For each biomolecular structure (e.g. PDB file), HL_viz.py (in Spectral-Eigvec folder) computes the Hodge Laplacians and generates the spectral eigenvectors via eigendecomposition.
```python
def faces(simplices):
    # This function outputs the faces of the input simplices. 
    
def n_faces(face_set, n):
    # Given a set of simplices as input, this function outputs the set of simplices of dimension n or of vertex set length n+1. 

def boundary_operator(face_set, i):
    # This function takes an input of set of simplices and output the boundary operator B_i.

# Input of structural data here
data = np.load("guanine.npz", allow_pickle=True)
for i in data["PRO"]:
    coords = i["pos"]

# Constructs the Vietoris-Rips complex at given cutoff distance e.g. 1.2.
rc = gd.RipsComplex(coords, max_edge_length=1.2)
simplex_tree = rc.create_simplex_tree(max_dimension=2)
val = simplex_tree.get_filtration()
simplices = set()
for v in val:
    simplices.add(tuple(v[0]))

# Computes the HL matrix by multiplying the boundary operators. 
laplacian = np.matmul(boundary_operator(simplices, 2).toarray(), np.transpose(boundary_operator(simplices, 2).toarray()))+np.matmul(np.transpose(boundary_operator(simplices, 1).toarray()), boundary_operator(simplices, 1).toarray())
eigval, eigvec = np.linalg.eigh(laplacian) # Performs the eigendecomposition.
print(eigval)
```

The magnitudes of eigenvectors can be stored in GEXF format via NetworkX which can then be viewed using GEPHI software. 
Examples are shown below.

![image](https://github.com/ExpectozJJ/Hodge-Theory/blob/main/images/guanine_eigenvec.png)

![image](https://github.com/ExpectozJJ/Hodge-Theory/blob/main/images/1axc_zero_eigenvec.png)

## Hodge Rank

The Hodge Decomposition codes are used for the following:

* DNA and Chromatin Folding
```
oglionucleosome_folding.m --> Computes the coboundary matrices, and the optimization to approximate gradient flow, harmonic flow and curl flows.
```
* Hi-C Data
```
HiC_oneDistance.m --> Computes the coboundary matrices, and the optimization to approximate gradient flow, harmonic flow and curl flows for the TAD Domains.
```
* Protein Folding
```
simulation_bending_protein.m --> Computes the coboundary matrices, and the optimization to approximate gradient flow, harmonic flow and curl flows for protein folding frames.
```

## Cite
If you use this code in your research, please cite our paper:

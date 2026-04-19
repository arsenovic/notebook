---
jupyter:
  jupytext:
    cell_metadata_filter: -all
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.18.1
  kernelspec:
    display_name: arsenovic-notebook (3.12.9)
    language: python
    name: python3
---

```python
# create  method to make a n-dimensional positive def symmetric amtrix 
import numpy as np
import kingdon as king 
import scipy 

def symmetric_matrix(n):
    A = np.random.rand(n, n)
    A = (A + A.T) 
    #A += n * np.eye(n)  # make it positive definite
    return A

def skew2bivector(A, alg):
    return alg.bivector(G[np.tril_indices(G.shape[0], k=-1)])  

n = 3
alg = king.Algebra(n)

Sigma = symmetric_matrix(n)
w, v  = np.linalg.eigh(Sigma)
G     = scipy.linalg.logm(v)
B     = skew2bivector(G,alg)
d     = alg.vector(w)
d+B
```

```python
scipy.linalg.logm(vecs).round()
```

```python


```

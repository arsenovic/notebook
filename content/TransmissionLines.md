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
from kingdon import Algebra
from kingdon.numerical import exp  as exp_
exp = lambda x: exp_(x,n=50) # precision control of exp

from kingdon.calculus import d
from kingdon.blademap import BladeMap
import numpy as np 
import torch
from timeit import default_timer

# e0 is dual origin, e4 = time , e123 is space
stap = Algebra(signature = [0,1,1,1,-1],start_index=0)
pga  = Algebra(3,0,1)
bm   = BladeMap(alg1=stap, alg2 =pga) 

S    = stap.blades
locals().update(S)
I    = stap.pseudoscalar([1])
 


 # bivector algebra
R = e3 * e4 - e1 * e3
X = -e2 * e4 + e1 * e2
G = e3 * e4 + e1 * e3
B = e2 * e4 + e1 * e2
N = e1 * e4
Q = e3 * e2
L = e1 * e2
A = e3 * e4


R**2, X**2, G**2, B**2, N**2, Q**2, L**2, A**2
```

```python
n= int(1e3)
r,x,g,b = 1,2,3,4
dr,dx,dg,db = r/n, x/n, g/n, b/n

lump = exp(r*R)*exp(x*X)*exp(g*G)*exp(b*B)  # a transmission line is not a lumped element
cell = exp(dr*R+dx*X +dg*G+ db*B)
cell**n , lump



```

```python

```

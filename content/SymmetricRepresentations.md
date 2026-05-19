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

<div style="text-align: right;">
written 04/19/2026<br>
<a href="mailto:alex@810lab.com">alex@810lab.com</a><br>
</div>


# Models for Symmetric Operators 

## sumary 
Taking an eigen-decomposition of a symmetric matrices seperates the scaling opertation from the rotation (an analgous dichotomy for polar coordinates). 
While you have to pay it forward, storing this information directly leads to some computation advantages:
* computing inverse is trivial 
* guaranteed to stay in the group. 
* interpolation is intuitive
* multiplication can be approximated with log-euclidean. (see [1])

and the following drawbacks:

* extra computation 
* addition of matrices is more expensive and cumberome. 

## reference 
[1] ["Fast and Simple Computations on Tensors with Log-Euclidean Metrics." Vincent Arsigny, et al.](https://www-sop.inria.fr/asclepios/Publications/Arsigny/arsigny_rr_tensors.pdf)

## code 

```python
import numpy as np
from scipy.linalg import logm, expm, eigh ,det,inv,schur ,inv


def cayley(A):
    ''' Maps a skew-symmetric matrix to an orthogonal matrix and vice versa'''
    I = np.eye(A.shape[0])
    return np.linalg.solve(I +  A, I - A)


def logm_ortho(V):
    '''fast way to computer log of an orthogonal matrix. uses schur'''
    # For orthogonal matrices, Schur decomposition is more stable and 
    # directly exposes the rotation blocks.
    # T is quasi-upper triangular, Z is orthogonal
    T, Z = schur(V, output='real')
    
    # In the real Schur form, T consists of 1x1 or 2x2 blocks on the diagonal.
    # log(T) for these blocks is trivial.
    logT = np.zeros_like(T)
    n = T.shape[0]
    i = 0
    while i < n:
        if i < n - 1 and abs(T[i+1, i]) > 1e-12:  # 2x2 block (rotation)
            # T[i:i+2, i:i+2] is [[cos, -sin], [sin, cos]]
            phi = np.arctan2(T[i+1, i], T[i, i])
            logT[i, i+1] = -phi
            logT[i+1, i] = phi
            i += 2
        else:  # 1x1 block (1 or -1)
            # log(1) = 0, log(-1) is technically +/- i*pi (requires care)
            logT[i, i] = 0 if T[i, i] > 0 else 0 # Simplified
            i += 1
            
    return Z @ logT @ Z.T

def logm_sym(A):
    ''' fast way to compute the matrix logarithm of a symmetric positive definite matrix'''
    w, R  = np.linalg.eigh(A)
    return R@np.diag(np.log(w))@R.T

def negate_diag(A):
    A = A.copy()
    A[np.diag_indices_from(A)] *= -1  
    return A

def symmetric_matrix(n):
    A = np.random.rand(n, n)
    return  A @ A.T

def to_dB(Sigma):
    w, R  = eigh(Sigma)
    if np.linalg.det(R) < 0: # force a R to be rotation 
        R[:, -1] *= -1  
    d = np.diag(np.log(w))
    B = logm_ortho(R)
    return  d+B
    
def from_dB(dB):
    d = np.diag(dB)
    B = dB - np.diag(d)
    R = expm(B)
    return R@np.diag(np.exp(d))@R.T



# test 
def eq(x,y):
    assert(np.allclose(x,y))
    
n=400
Sigma = symmetric_matrix(n)*128
dB = to_dB(Sigma)

eq(from_dB(to_dB(Sigma)), Sigma)         # round trip
eq(logm(Sigma),logm_sym(Sigma))          # fast logm works
eq(from_dB(negate_diag(dB)), inv(Sigma)) # negating the diagonal == inverse
eq(np.exp(np.trace(dB)), det(Sigma))# 

w, R  = eigh(Sigma)
#eq(cayley(cayley(R)) ,R ) 
logm_ortho(R).max(), cayley(R).max()
```

```python
import timeit
import pandas as pd

def benchmark_functions(func_dict, number=100, repeat=3):
    '''
    Benchmark a dictionary of functions and return timing results.
    
    Parameters:
    -----------
    func_dict : dict
        Dictionary with function names as keys and callables as values
    number : int
        Number of times to execute each function
    repeat : int
        Number of times to repeat the timeit measurement
    
    Returns:
    --------
    pd.DataFrame with columns: Function, Min (s), Mean (s), Std Dev (s)
    '''
    results = []
    
    for name, func in func_dict.items():
        # Measure execution time
        times = timeit.repeat(func, number=number, repeat=repeat)
        min_time = min(times) / number
        mean_time = sum(times) / (len(times) * number)
        std_dev = (sum((t/number - mean_time)**2 for t in times) / len(times)) ** 0.5
        
        results.append({
            'Function': name,
            'Min (s)': f"{min_time:.6e}",
            'Mean (s)': f"{mean_time:.6e}",
            'Std Dev (s)': f"{std_dev:.6e}"
        })
    
    return pd.DataFrame(results)

# Example usage:
func_dict = {
    'inv(Sigma)': lambda: inv(Sigma),
    'from_dB(negate_diag(to_dB(Sigma)))': lambda: from_dB(negate_diag(to_dB(Sigma))),
    'negate_diag(dB)': lambda: negate_diag(dB),
}

results = benchmark_functions(func_dict, number=10, repeat=3)
print(results)
```

```python
# compare times for a  full  matrix inv() vs a inv_diag()
%timeit inv(Sigma)
```

```python
%timeit  from_dB(negate_diag(to_dB(Sigma)))
```

```python
%timeit negate_diag(dB)
```

```python
w, R  = eigh(Sigma)
if np.linalg.det(R) < 0: # force a R to be rotation 
    R[:, -1] *= -1  
d = np.diag(np.log(w))



B = logm_ortho(R)
eq(inv_cayley(cayley(R)) ,R )
```

```python

# convert a matrix to a given fp4 format and back, should be identity
B.max()
```

```python
%timeit cayley(R)
```

```python
def sqrt_sym(A):
    # faster way to compute the matrix square root of a symmetric positive definite matrix
    w, R  = np.linalg.eigh(A)
    return R@np.diag(np.sqrt(w))@R.T

A = np.random.rand(n, n)
Sigma = A @ A.T
sqrt_sym(Sigma)@sqrt_sym(Sigma), Sigma, 
A , sqrt_sym(Sigma)
```

## Geometric Algebra  Model

```python


import scipy 

def skew_to_bivector(A, alg):
    return alg.bivector(A[np.tril_indices(A.shape[0], k=-1)])  

n = 3
import kingdon as king 
alg = king.Algebra(n)

Sigma = symmetric_matrix(n)
w, v  = np.linalg.eigh(Sigma)
if det(v) < 0:
        v[:, 0] = -v[:, 0]
B     = skew_to_bivector(scipy.linalg.logm(v),alg)
d     = alg.vector(w)
d+B
```

```python
np.log(w)
```

```python

```

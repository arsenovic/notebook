---
jupyter:
  jupytext:
    default_lexer: ipython3
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.18.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

<!-- #region -->
# In Search of  Nilpotents
alex arsenovic 6/19/25
## Intro 

This writeup is a work in progress toward understanding nilpotent matrices in Geometric Algebra . 


Symmetric matrices are important, but the GA spinor-model for symmetric matrices is awkward.  the current spinor model  requires 2n-dimensions, and implements anistropic scaling with  hyperbolic rotations of   null vectors,   see Doran's balanced algebra. In addition, the cholesky decomposition has come up a lot, due to its numerical advantages. 

### Cholesky 
So, we began looking for a GA model for Cholesky, which states that a symmetric matrix $S$ can be written as the product of two triangular matrices, $L$ and $L^T$ .   

$$ S = LL^T $$

the $L$ can also be further decomposed into,
$$ L = D e^N $$

(aka LDL), where $D$ is diagonal and $N$ is nilpotent of degree $k$, meaning $N^k = 0$ while $N^j \ne 0$ for $1 \le j < k$. 
I think  this can be interepreted as  decomposing $S$ into shears and strains. (see CA2GC 3-6)

### Nilpotent
The nilpotent property is interesting, 
* the solution to a [matrix differential equation](https://en.wikipedia.org/wiki/Matrix_differential_equation) is $e^{tN}$ for  polynomial (kinematic) state model
*  can be used to take derivatives of polynomials (see autodiff).
  
The nilpotent  property makes the matrix exponential function expression terminate after $k$-terms 

$$
\exp N \;=\; I+N+\tfrac12 N^{2}+\tfrac16 N^{3},
$$


i am not sure how to best model this either. what multivector could have this property? 
Lets ignore the strains ($D$) for now, since we have failed at this many times over, and look for the ability to represent nilpotents $N$ in $G_n$.
<!-- #endregion -->

### Shift matrix 
Any nilpotent matrix  is similar to a a block-diagonal of shift matrices [wiki](https://en.wikipedia.org/wiki/Shift_matrix) . This seems like the ticket. A shift matrix is  a matrix with ones along the super diagonal. 

$$N_5 = \begin{pmatrix}
0 & 1 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 1 \\
0 & 0 & 0 & 0 & 0
\end{pmatrix}$$
 as you square it, the super diagonal  moves away from the diagonal, and eventually you get $N_5^5=0$ .


#### Functional form 
its easy to write a shift matrix as a geometric  function, being a sum of shears. 
$$\sum_{i=1}^{k} \alpha_i\, (x\cdot e_i)\, e_{i+1}$$

This is very  similar to the diagonal matrix , being a sum of strains,

$$\sum_{i=1}^{k} \alpha_i\, (x\cdot e_i)\, e_i$$

The inner product can be expanded, to give an alternative interpretation for both expressions


### Circular Shift 
 A closely related matrix is a circular shift.  this can obviously be written as a product of reflections using Cartan's algo. 
 $$C_5 = \begin{pmatrix}
0 & 1 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 1 \\
1 & 0 & 0 & 0 & 0
\end{pmatrix}$$

but, not clear to reduce it to the shift matrix. i thought making on dimension null ($e_i^2=0$), might work but it didnt seem to.


### Sum of positive and negative bivectors?

would this sort of diachotomy be useful? 
$$N= \frac{1}{2}(N_5^{+}+N_5^{-} )$$
where

$$N_5^{+} = \begin{pmatrix}
0 & 1 & 0 & 0 & 0 \\
1 & 0 & 1 & 0 & 0 \\
0 & 1 & 0 & 1 & 0 \\
0 & 0 & 1 & 0 & 1 \\
0 & 0 & 0 & 1 & 0
\end{pmatrix}$$
$$N_5^{-} = \begin{pmatrix}
0 & 1 & 0 & 0 & 0 \\
-1 & 0 & 1 & 0 & 0 \\
0 & -1 & 0 & 1 & 0 \\
0 & 0 & -1 & 0 & 1 \\
0 & 0 & 0 & -1 & 0
\end{pmatrix}$$


## End
Thats all we got for now


## Code
numerical tests of functions and how they effect a set of basis vectors

```python
def basis_vectors_lst(alg):
    return [alg.blades.blades[k] for k in alg.blades.blades.keys() if len(k)==2]

def map_check(func,alg):
    a = basis_vectors_lst(alg)
    for k in range(len(a)):
        print(a[k], "->", func(a[k]))

def func_2_mat(func,alg):
    A = basis_vectors_lst(alg)
    B = [func(a) for a in A]
    M = np.array([(a|b).e for a in A for b in B]).reshape(len(A), len(A)) 
    return M

def shift(v):
    return sum( a[k-1]*(a[k] | v)   for k in range(1, n) )

def circ_shift(x):
    rs=[ (a[k]+a[(k+1)%n]).normalized()  for k in range(n)]
    R = rs[0]
    for r in rs: R*=r
    #R = math.prod(rs)#?
    return R>>x

#####
import numpy as np
import kingdon as king 
import math 

n   = 4
alg = king.Algebra(n,0,0)
#alg = king.Algebra(n-1,0,1)
a   = basis_vectors_lst(alg)

f = shift
f = circ_shift
map_check(f, alg), func_2_mat(f,alg)
```

### Inital Brute-force


We did some chatgpt'ing and we know that a rotation in null bivector implements a shear. so we guessed to use G(N,0,1) and brute force found some nilpotents as show below.  the patern and meaning is to be found.

```python
from kingdon import Algebra
alg = Algebra(3, 0,1)
locals().update(alg.blades)

X = e2*(1 + e1*(1 + e0))
[print(f'X^{n} : {(X)**n}') for n in range(1,5)];
```

### G(3,0,1)

```python
from kingdon import Algebra
alg = Algebra(3, 0,1)
locals().update(alg.blades)

import math 
def exp_nil(X):
    return sum([X**k/math.factorial(k) for k in range(X.algebra.d)])

def log_nil(B):
    K = B-1
    return sum([(-1)**(k+1) * K**k/k for k in range(1,B.algebra.d)])

N = e2*(1 + e1*(1 + e0))
assert((log_nil(exp_nil(N))-N).e ==0)

[print(f'N^{n} : {(N)**n}') for n in range(1,5)];
```

```python
Ns = [m*(1 + n*(1 + e0)) for m,n in [(e2,e1),(e2,e3),(e1,e3)]]
[(k**2,k**3) for k in Ns]
```

```python
Ns[0].cp(Ns[1]), Ns[2]
```

```python
Ns
```

```python
N1,N2,N3 = Ns

R1= exp_nil(N1)

R1*R1
```

### G(4,0,1)

```python
from kingdon import Algebra
alg = Algebra(4, 0,1)
locals().update(alg.blades)

Z = e2*(1 + e1*(1 + e0))
Y = e4*(1 + e3*(1 + e0))
X = Z + e0*Y
[print(f'X^{n} : {(X)**n}') for n in range(1,5)];
```

```python
N
```

```python

```

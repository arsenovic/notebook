---
jupyter:
  jupytext:
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

## Linear Geometric Functions

Some classes for interepreting linear operators and converting to geometric functions. 

```python code_folding=[11, 28, 43, 66, 71, 90]

from functools import reduce,wraps
import warnings
from clifford import op 
from clifford.tools import orthoMat2Versor,log_rotor,mat2Frame,frame2Mat#,func2Mat
#from clifford.invariant_decomposition import bivector_split,log
import numpy as np
from scipy.constants import e, pi 
from typing import Union, Optional, List, Tuple
cayley =  lambda x:(1-x)/(1+x)

def func2Mat(f,I):
    '''
    Convert a function to a matrix by acting on standard basis

    Parameters
    ---------------
    f : function
        function that maps vectors to vectors
    I : MultiVector
        psuedoscalar of basis

    See Also
    ---------
    frame2Mat
    '''
    A = I.basis()
    B = [f(a) for a in A]
    return frame2Mat(B=B, A=A,I=I)

def outermorphic(f):
    '''
    convert a geometric function `f` which operates on vectors into a 
    geometric function which operates on multivectors, using outermorphism.
    
    useful as decorator
    '''
    def F(M): # M is a multivector
        N = []
        for blade in M.blades_list:            # decompose MV into blades
            factors,scale  = blade.factorise() # decompose each blade into factors 
            y = [f(x) for x in factors]        # f(x) each factor 
            Y = reduce(op,y)*scale             # compile factors back into blades
            N.append(Y)                        # compile blades back into MV
        return sum(N)
    return F 

def symmetric_2_skewup1(M):
    '''convert a symmetric matrix of rank `n` to a skewmetric matrix of rank `n+1`'''
    rank = M.shape[0]
    out= np.zeros([rank+1,rank+1])
    out[0,1:]= np.diagonal(M)
    out[1:,1:]=M
    out = np.triu(out,1)
    return out-out.T

class MatrixOperator(object):
    def __init__(self,M,I=None):
        '''
        An abstract base class inhereted by matrix operators
        '''
        if I is None: # make a euclidean N-space
            layout,blades = Cl(M.shape[0])
            I = layout.pseudoScalar
            
        self.I = I 
        self.M = M
        self.validate_M()   
    
    @property
    def rank(self):
        return self.M.shape[0]
    
    def validate_M(self):
        pass
      
class Linear(MatrixOperator):
    '''
    A linear operator, as defined by a square matrix. 
    '''
    @classmethod
    def from_rand(cls, n):
        return cls(M = np.random.randint(-100,100,n**2).reshape(n,n))
        
    @property
    def Me(self):
        return (self.M + self.M.T)/2.
    
    @property
    def Mo(self):
        return (self.M - self.M.T)/2.
    
    @property
    def symmetric(self):
        return Symmetric(I= self.I ,M = self.Me)
    
    @property
    def skewmetric(self):
        return Skewmetric(I= self.I ,M = self.Mo)
    
    def as_f(self):
        raise NotImplementedError('works if self.normal?') 
        g = self.symmetric.as_f()
        h = self.skewmetric.as_f()
        return lambda x: g(x) + h(x)
    
class Skewmetric(MatrixOperator):
    '''
    A skew-symmetric (aka skewmetric)  operator. 
    '''
    def validate_M(self):
        if not np.allclose(self.M,-self.M.T):
            warnings.warn('M is not skewmetric')
        
    def as_bivector(self):
        '''
        convert a skewmetric matrix to its curling bivector
        '''
        B = mat2Frame(self.M,I=self.I)[0]
        A = mat2Frame(np.eye(self.rank),I=self.I)[0]
        F = sum([a*b for a,b in zip(A,B)])/2
        return F

    def as_f(self):
        '''
        convert a skewmetric (antisymmetric) matrix into a
        geometric function `f`
        '''
        F   = self.as_bivector()
        f_x = lambda x: x|F
        f   = outermorphic(f_x)
        return f
    
class Symmetric(MatrixOperator):
    '''
    A symmetric operator
    '''
    def validate_M(self):
        if not np.allclose(self.M,self.M.T):
            warnings.warn('M is not symmetric')
            
    def as_eigenframe(self):
        '''
        Return the frame (list of vectors) for the eigen vectors, where the 
        eigenvectors are scaled by their corresponding values
        '''
        w,v = np.linalg.eig(self.M)
        return mat2Frame(v@np.diag(w),I=self.I)[0]
    
    def as_vector_and_versor(self):
        w,v = np.linalg.eig(self.M)
        d   = sum([a*k for a,k in zip(w,self.I.basis())])
        R,rs  = orthoMat2Versor(v,I= self.I)
        return d,R
        
    def as_f(self):
        '''
        convert a symmetric matrix into one possible geometric function `f`
        '''
        #V = self.symmetric_mat_2_eigenframe(M)
        #f_x = lambda x: sum([(x|k)/k*abs(k) for k in V])
        w,v = np.linalg.eig(self.M)
        V   = mat2Frame(v,I=self.I)[0]
        f_x = lambda x: sum([(x|k)/k*a for k,a in zip(V,w)])
        f   = outermorphic(f_x)
        return f
    
    def as_skewup1(self):
        '''
        convert a symmetric operator into a skewmetric operator in 
        a space one dimension higher. 
        '''
        M = symmetric_2_skewup1(self.M)
        layout,blades = Cl(M.shape[0], firstIdx=0)
        I = layout.pseudoScalar
        
        return Skewmetric(M=M,I=I)


```

```python
from clifford import Cl

n = 3
M = np.random.randint(-100,100,n**2).reshape(n,n)

layout,blades = Cl(n)
I      = layout.pseudoScalar
linear = Linear(I=I, M=M)
#linear.skewmetric.as_bivector()

```

```python
from clifford import Cl

n = 2
M = np.random.randint(-100,100,n**2).reshape(n,n)

layout,blades = Cl(n)
I      = layout.pseudoScalar
linear = Linear(I=I, M=M)
linear
```

```python
from clifford import Cl

n = 2
M = np.random.randint(-100,100,n**2).reshape(n,n)

layout,blades = Cl(n)
I      = layout.pseudoScalar
linear = Linear(I=I, M=M)


assert(np.allclose(func2Mat(linear.symmetric.as_f(),  I=I)[0], linear.Me))
assert(np.allclose(func2Mat(linear.skewmetric.as_f(), I=I)[0], linear.Mo))
a,b = linear.symmetric.as_eigenframe()
assert(a|b==0)
```

```python
linear = Linear.from_rand(2)

diag = np.diag([1,2])
R = linear.I.layout.randomRotor()
Rmat,dum = func2Mat(lambda x: R*x/R, I = linear.I )

M = Rmat@diag@Rmat.T
linear = Linear(M)

linear.M

```

```python
linear.symmetric.M
```

```python
a,b = linear.symmetric.as_eigenframe()
abs(a), abs(b)
```

```python
f = linear.symmetric.as_skewup1().as_f()
F = linear.symmetric.as_skewup1().as_bivector()
F
```

```python
locals().update(linear.symmetric.as_skewup1().I.layout.blades)

f(e1)
```

```python
f(F)/F
```

```python
f = linear.skewmetric.as_f()
F = linear.skewmetric.as_bivector()
F





```

```python
linear.symmetric.M
```

```python


```

```python
f   = linear.skewmetric.as_f()
F   = linear.skewmetric.as_bivector()
Fks = bivector_split(F)
[f(Fk)/Fk for Fk in Fks] # should all be scalar. def of skewsymmetry and invariant eigenbivector
```

```python
Fks
```

```python
linear.skewmetric.M
```

```python
linear.skewmetric.as_bivector()
```

```python
n=3
M       = np.random.rand(n,n)
M       = .5*(M+M.T)
symmetric_2_skewup1(M),M
```

```python
linear.Mo+linear.Me,linear.M
```

```python code_folding=[]
from clifford import Cl
n=3

lay,b = Cl(sig=[1]+list(np.ones(n)),firstIdx=0)
locals().update(b)
M       = np.random.rand(n,n)
M       = .5*(M+M.T)
lo      = Symmetric(I=lay.pseudoScalar*e0, M=M)
hi      = Skewmetric(I=lay.pseudoScalar, M=symmetric_2_skewup1(lo.M))

d,G = lo.as_vector_and_versor()
#(d*e0),log_rotor(G)
d,G
```

```python
import scipy 
M = np.random.rand(n,n)
l = Linear(M=M,I=I)
scipy.linalg.polar(np.random.rand(n,n)), l.symmetric.M


```

```python
,F = hi.as_bivector()
G = log_rotor(cayley(F))

```

```python
f = hi.as_f() 
Ks = [k*e0 for k in lo.as_eigenframe()]

[f(K)/K for K in Ks]

```

```python
Cl(sig=[-1]+list(ones(2)),firstIdx=0)
```

```python
## commutator patterns
def symmetric_M(n):
    M     = np.random.randint(0,20,n**2).reshape(n,n)
    M = (M+M.T)/2.
    return M
def skewmetric_M(n):
    M     = np.random.randint(0,20,n**2).reshape(n,n)
    M = (M-M.T)/2.
    return M

M1 = symmetric_M(3)    
M2 = symmetric_M(3)    
N1 = skewmetric_M(3)
N2 = skewmetric_M(3)

com = lambda x,y: x@y-y@x

```

```python

```

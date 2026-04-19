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
## Whitening  

### Summary
Here, we model a gaussian distribution as an Affine operator. This is done by  using homogeneous coodinates which combines the mean and covariance into a single matrix. This leads to a compact representation  and intuitive  interpretation of the mahalonobis distance, pdf, etc. 

This concept was done in the 90s, see  "A distance between multivariate normal distributions based in an embedding into the siegel group". by  Miquel Calvo, Josep M. Oller.

### Context
This idea came from the wiki page on statistical [whitening](https://en.wikipedia.org/wiki/Whitening_transformation): 

"*A whitening transformation or sphering transformation is a linear transformation that transforms a vector of random variables with a known covariance matrix into a set of new variables whose covariance is the identity matrix, meaning that they are uncorrelated and each have variance 1.[1] The transformation is called "whitening" because it changes the input vector into a white noise vector.*"  



### Setup
A multivariate Gaussian distribution in '$k$' dimensions is given by 

$$pdf(x) = \frac{1}{\sqrt{(2\pi)^k |\Sigma|}}
e^{(x-\mu)^T\Sigma^{-1}(x-\mu)}$$

the term in the exponent is known as the 'mahalonobis distance' which is basically a distance in units of standard deviations.

$$d^2_m(x) =(x-\mu)^T\Sigma^{-1}(x-\mu) $$
It can be re-written as a norm  by interpretting the covariance $\Sigma$ as an operator, and splitting the operation into two halves, loosely speaking. There are several canonical ways to decompose a matrix (SVD, cholesky, etc).   One way is 
$$\Sigma=AA^T$$
This decomposition will be Cholesky if $A$ is  triangular, but we dont need to make that choice yet. 
#### Side Notes

* Notice that in 1-D, variance is written $\sigma^2$, while in n-D its just $\Sigma$. However, the units of covariance are that of $\Sigma^2$, so i would call this a math nomenclature bug. The fact that it should be $\Sigma^2$ actually suggests the decomposition,  shown above also. 

* For any matrix $A$, the product $AA^T$ is symmetric. For *any* matrix  $A$. crazy.

In either case, the inverse of $\Sigma$  (known as the precision matrix) is  given in terms of $A$ as  
$$\Sigma^{-1}=(A^T)^{-1}A^{-1} = (A^{-1})^{T}A^{-1}. $$

The maha distance expressed in terms of this half operation $A$ is
<!-- #endregion -->

<!-- #region -->
\begin{align}
d^2_m(x) &=(x-\mu)^T\Sigma^{-1}(x-\mu) \\
&= (x-\mu)^T(A^{-1})^{T}A^{-1}(x-\mu)\\
&=|A^{-1}(x-\mu)|^2 
\end{align}


Which can be interpreted as the norm of a vector after the  inverse of some affine transform, call it $g$. 
$$  g(x) = Ax+\mu, \qquad g^{-1}(x) = A^{-1}(x-\mu),$$
$g$  is a combination of a linear transform $A$, followed by a translation of $\mu$. This suggests the interpretation of $g$ as a 'whitening' operator which transforms between a *standard* normal  and the distribution of interest. In other words,  it is an operator model for a distribution. This is analagous to the way  complex numbers are used in LTI system analysis, (as operators on a reference sinusoid) . 

### Matrix Model 
 The operation (linear+translation) can be implemented with a single affine transformation matrix $G$ using [homogenous coordinates](https://en.wikipedia.org/wiki/Affine_transformation#Augmented_matrix),  
$$
g(x) = \downarrow G \mathbf{x} = \downarrow \underbrace{ \left[ \begin{array}{ccc|c} & A & & \mathbf{\mu} \\ 0 & \cdots & 0 & 1 \end{array} \right]}_{G}
\underbrace{\begin{bmatrix} x \\ 1 \end{bmatrix}}_{\mathbf{x}}
$$
and 
$$
g^{-1}(x) = \downarrow G^{-1}\mathbf{x}
$$
 
Where $\downarrow$ implies downprojection from affine space/homogeneous coordinates. In terms of the affine matrix $G$, The maha distance is,
$$d_m^2 =|g^{-1}(x)|^2= 1-|G^{-1}\mathbf{x}|^2$$
(the '$1-$' is just a little hack to make it work with the homogeneus coordinates which avoids the $\downarrow$, there are other ways to do this too). The formula for the pdf in terms of homogenous coordinates becomes, 

$$\text{pdf}(\mathbf{x}) = \frac{1}{(2\pi)^{k/2} |G|}
e^{1-|G^{-1}\mathbf{x})|^2}$$
(note $|G| = \sqrt{|\Sigma|}$), or put another way 
$$\text{pdf}(\mathbf{x}) = \frac{e}{(2\pi)^{k/2}}|G^{-1}|
e^{-|G^{-1}\mathbf{x}|^2}$$

So in some sense, we can say that $G$ *generates*  the pdf. kindof.
<!-- #endregion -->

### Propagation

Normally, if we want to propagate a distribution through a linear operator,  you have two rules to implement; one for the covariance and another for the mean. 
$$ \Sigma' = F\Sigma F^T, \qquad \mu'  = F\mu$$
but since
$$ F\Sigma F^T = FAA^TF^T  $$
The propagation of $A$ is the same as $\mu$. So in the affine representaion its  just one rule ( and one matrix multiply)
$$ G' = FG$$

Now, while this is nice,  it presents a practical problem if you want to get the covariance matrix out, or add two covariance matrices ([see here](https://mathoverflow.net/questions/364442/how-do-you-find-the-cholesky-decomposition-of-the-sum-of-two-positive-definite-m)). In this case  you have to form  $AA^T$.


Lets implement  this. Note there are some annoying book-keeping depending on how you decompose $\Sigma \rightarrow AA^T$ or $A^TA$, and so on.

```python
import numpy as np 
from numpy import linalg as la
from numpy import e 
import math

def rand_vec(n): # random vector 
    return np.random.random(n)

def rand_mat(n): # random matrix 
    return np.random.random((n,n))

def rand_spd(n): # random symmetric positive definite matrix 
    A = np.random.random((n,n))
    return A@A.T

def affinize(S=None,x=None):
    ''' make affine matrix/homogenous vector from linear operator and/or translation vector '''
    if S is None and x is not None:
        return np.append(x,1) # return a homogenous vector 
    else: # return a matrix 
        n = len(S)
        A = np.eye(n + 1)       
        A[:n, :n] = S
        if x is not None:
            A[:n, n]  = x.flatten()
    return A 

def afg(Sigma,mu=None):
    ''' create an affine matrix for a gaussian distribution'''
    L = la.cholesky(Sigma,upper=True)
    return affinize(S=L, x=mu)
    
################## Tests #############################
n=3

def test_homogenous():
    # homogenous matrix tests 
    S  = rand_spd(n)
    mu = rand_vec(n) 
    x  = rand_vec(n) 
    
    X  = affinize(x=x)
    A  = affinize(S,mu)
    Ai = la.inv(A)
    Si = la.inv(S)
    
    assert(np.allclose(( A@X)[:n],(S@x+mu))) # sanity checks
    assert(np.allclose((Ai@X)[:n], Si@x-Si@mu))
test_homogenous()

def test_affine_gaussian():
    F     = rand_mat(n)
    x     = rand_vec(n)
    Sigma = rand_spd(n)
    iSigma = la.inv(Sigma)
    mu    = rand_vec(n)
    
    X     = affinize(x=x)
    G     = afg(Sigma=Sigma, mu=mu)
    A     = G[:n,:n]
    
    iL    = la.cholesky(iSigma,upper=True)
    iG    = affinize(iL, -iL@mu)
    #why not iG    = la.inv(G)#afg(Sigma=la.inv(Sigma), mu=-mu)#la.inv(G)
    #np.allclose((G.T@G)[:n,:n] , Sigma)
    
    assert(np.allclose(A.T@A, (G.T@G)[:n,:n], Sigma))
    assert(np.allclose((x-mu).T@iSigma@(x-mu), (iG@X)@(iG@X)-1) )# maha tests 
    assert(np.allclose( la.det(Sigma), la.det(G)**2 , math.prod(np.diag(G))**2)) # det tests
    FGt = affinize(F)@G.T
    assert(np.allclose( (FGt@FGt.T)[:n,:n], F@Sigma@F.T)) # function tests 

test_affine_gaussian()
```

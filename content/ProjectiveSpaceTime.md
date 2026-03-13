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
    display_name: arsenovic-notebook
    language: python
    name: python3
---

<!-- #region -->
Here we re-interpret a plane wave in terms of the projective geometric algebra.

Consider  the form of a EM plane wave given in Hestenes' Spacetime Algebra 
$$ F_x \equiv  F(x)= e^{I(k\cdot x)}F_0$$

where;
* $x$ is a location in spacetime
* $k$ is the propagation vector
* $I$ is psuedoscalar. 

The generator re-expressed in a couple ways. 

$$  I(k\cdot x) =K\wedge x =   I(k\cdot x) =  k\wedge X  $$

Through the lense of PGA, we can interepret  $X=Ix$ as a point in spacetime, and $k$ a  plane. Then k\wedge X is the signed distance between a point and plane, which parameterizes the duality rotation. 

From this perspective time-harmonic analysis requires  $\omega \equiv k\cdot \gamma_0 \gamma_0^{-1}   $ to be a constant. This produces a projective geometry on the $\omega$ plane.  

[Previously](PlaneWaves.html) we considered interpretting the generator as,  



This equation is interpreted as a duality rotation of a fixed bivector $F_0$, parameterized by $K$, with $x$ being  an independent variable.  The vector derivative of this is ,
$$ \nabla e^{K\wedge x} = K e^{K\wedge x}   $$

 
<!-- #endregion -->

```python

from kingdon import Algebra
from kingdon.numerical import exp , is_close
from kingdon.calculus import d,curl,div


sta = Algebra(signature=[-1,-1,-1,1],start_index=1)
locals().update(sta.blades)
I = sta.pseudoscalar([1])
D = lambda f, **kw: d(sta, f, **kw) # vector derivative in sta

def round(mv, tol=1e-6):
    return mv.filter(lambda c: abs(c) > tol)

K  = sta.random_trivector() 
k  = I*K
x  = sta.random_vector()
B  = sta.random_bivector()
a0 = K|B



assert is_close(I*(k|x) , x^K)
f = lambda x:  exp(x^K)
assert is_close(D(f)(x), K*f(x)) 
#D(D(f))(x), -K*K*f(x) # why this not work?
assert is_close(D(lambda x: exp(I*(k|x)))(x),-I*k*exp(I*(k|x))) # same as above

Kf = lambda x: K*exp(x^K)
assert is_close(round(D(Kf)(x)), -K**2*f(x)) # why doe we need a round here?
```

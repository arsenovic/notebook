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

Here we attempt to write down a plane wave in terms of the vector potential with the least number of free parameters. 

Consider  the form of a EM plane wave given in Hestenes' Spacetime Algebra 

$$ F_x \equiv  F(x)= e^{I(k\cdot x)}F_0=  e^{K\wedge x}F_0 $$

where;
* $x$ is a location in spacetime
* $k$ is the propagation vector
* $K=Ik$ is a trivector dual to $k$.
* $I$ is psuedoscalar. 

This equation is interpreted as a duality rotation of a fixed bivector $F_0$, parameterized by $K$, with $x$ being  an independent variable.  The vector derivative of this is ,
$$ \nabla e^{K\wedge x} = K e^{K\wedge x}   $$

The electromagnetic bivector can be expressed as the curl of some vector field  $A_x$  
$$\nabla \wedge A_x = F_x$$

Lets   consider a vector derivative-relationship
$$\nabla  A_x = F_x$$

Then, if $F_x$'s x-dependence is a duality rotation, then so must $A_x$'s ,
$$A_x = e^{K\wedge x}A_0$$

which implies
$$\nabla A_x = F_x = e^{K\wedge x}KA_0.$$

(Moving the $K$ through the duality rotation at most picks up a sign change that we ignore). So, 
$$F_0 = KA_0  $$

For this $F_0$ to be a bivecotr , $K \wedge A_0 =0$. 
One way to solve this for $A_0$ in terms of $K$ is to pick a random bivector $B$, and contract it with $K$
$$A_0 = K\cdot B = \langle KB \rangle_1$$
then 
$$ F_0= K (K\cdot B)$$
Leading to the equation 
$$ F_x = e^{K\wedge x}K (K\cdot B)$$

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

```python
f = lambda x:  exp(x^K)

```

```python
Kf = lambda x: K*f(x)
round(D(Kf)(x)), -K**2*f(x)
```

```python
div(sta,div(sta, f))(x),curl(sta,curl(sta, f))(x)


```

```python



a = lambda x:  exp(x^K)*a0
F=d(sta,a)

a0, K, a0^K, K**2,round(F(x)), K|B


```

```python
sta.random_trivector()*sta.random_vector()
```

```python
d(sta,lambda x: exp(I*(k|x)))(x),-I*k*exp(I*(k|x))
```

```python
from kingdon import Algebra
from kingdon.numerical import exp 
from kingdon.calculus import d,curl
sta = Algebra(signature=[-1,-1,-1,1],start_index=1)
locals().update(sta.blades)
F = sta.random_bivector()
I = sta.pseudoscalar([1])
import numpy as np 

norm = lambda x: x/np.sqrt(abs(((x)**2).e))

def decomp(F, p,tol=1e-6):
    if abs((p**2).e) < tol:  # we have a null vector. use idempotents.
        w,k = decomp(p, e4)
        w,k = norm(w), norm(k)
        ep, em =(w+k)/2, (w-k)/2
        epm, emp  = -ep*em, -em*ep
        return F*epm, F*emp
    else:
        return 1/2*(F + p*F/p), 1/2*(F - p*F/p)

x = sta.random_vector()
#P = sta.random_trivector()
P = e123+e124

a = lambda x: exp(x^P)*e1




a(x)*e4, curl(sta,a)(x)
```

```python
exp(x^P)
```

```python
d(sta, lambda x: exp(x^P))(x), P*exp(x^P)
P*sta.random_multivector().odd

```

```python
P|e4
```

```python
W,K = decomp(P,e4)
((W/I)^(K/I))*I

```

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
# Plane Waves 

## summary 
An equation is reverse engineered from the STA EM plane wave  equation which follows conjugation. 


## Idea

Consider  the form of a EM plane wave given in Hestenes' Spacetime Algebra 

$$ F(x)= F_x = e^{I(k\cdot x)}F_0 $$

where;
* $x$ is a location in spacetime
* $k$ is the propagation vector
* $I$ is psuedoscalar. 

This equation is interpreted as a duality rotation of a fixed bivector $F_0$, parameterized by $k$, with $x$ being  an independent variable. It can also be written as 

$$ F_x = e^{(Ik)\wedge x}F_0 =  e^{K\wedge x}F_0$$

where $K$ is a trivector dual to $k$. ( We prefer trivectors  for propagation constants since the units match the grade ( rad/m rad/s), but keeping with dual vectors works just as good.) This function solves the equation, 
$$ \nabla F_x= KF_x. $$


<!-- #endregion -->

<!-- #region -->

So far we have changed nothing from STA. If we re-express $F_x$ in terms of the geometric product, it suggests to us a formula which employs conjugation, such as  

$$ e^{K\wedge x}F_0 = e^{Kx-xK}F_0 \stackrel{?}{=} e^{Kx}F_0e^{-xK}. $$

(we ignore all factors of 1/2 untill we need to precise answer).  For the last  equality to hold, $F_0$ must commute with the bivector $K\cdot x$. However, since $x$ is an independent variable, this cannot be guaranteed without more assumptions. So,  there are two possibilities 
1) this reverse engineered formula is useless
2) The plane-wave forumula is making implicit assumptions which justify the commutation assumption. 

For now, lets proceed and  see what such a forumla provides. Expanding the product 

$$xK = x\cdot K + x \wedge K  = B+\alpha I  $$

Since the psuedoscalar commutes with $F_0$, we can then write, 


$$  e^{ K\wedge x} e^{K\cdot x} F_0e^{-K \cdot x} = e^{ \alpha I }e^{B}F_0e^{-B}   $$

This is interepreted as a STA position-dependent lorentz rotation and a duality rotation. Along certain axis, namely $x\cdot K = 0$, the lorentz transforms disappears and this reduces to the normal plane wave.

### Derivative

<!-- #endregion -->


```python

```

## extensions

* trivector pregenerator becomes trivector+vector to model loss $K\rightarrow K+k$
* extend  $F$ to be a arbitrary multivector $M$
* trivector/vector pregenerator to also be a function of position, $K\rightarrow K_x$ 

$$ M_x = e^{K_x x}M_0e^{-xK_x}. $$

 ## some Details 
 * order of $xK$ vs $Kx$,  and + vs -
 * re-work in terms fo vector $k$ 




## Field Visualization 

```python
from kingdon import Algebra
from kingdon.numerical import exp
from kingdon.calculus import d
from kingdon.blademap import BladeMap
import numpy as np 
from timeit import default_timer
c = [0, 1810039, 14245634, 7696563, 15149450, 6727198, 15117058, 10909213, 6710886]  # colors 

sta = Algebra(1,3,start_index=0)
D = sta.blades
DI  = sta.pseudoscalar([1])


pga = Algebra(3,0,1,start_index=0)
P = pga.blades
bm = BladeMap(alg1=sta, alg2 =pga) 
up = lambda x: (x+P.e0).dual()




def graph_func():
    t = default_timer() / 100000
    #print(t)

    N=2
    x1, x2 = np.mgrid[-1:1:N*1j, -1:1:N*1j]
    x = x1*D.e1+x2*D.e2+t*D.e0
    xs= x.flatten()

    F0 = D.e01 + D.e02.dual # E1+H2
    K = D.e123
    F = lambda x: exp(K^x)*F0

    ## convert to cga to display 

    Fs = [F(x) for x in xs]
    E,H = zip(*[(F.proj(D.e0), F.proj(D.e0.dual).dual()) for F in Fs])
    E= [k.normalize() for k in E]
    H= [k.normalize() for k in H]
    E = list(map(bm,E))
    H = list(map(bm,H))
    Xs = list(map(up, (map(bm, xs))))
    vmag = .2
    ve = [(x,.5*exp(vmag*k)>>x) for x,k in zip(Xs,E)]
    vh = [(x,.5*exp(vmag*k)>>x) for x,k in zip(Xs,H)]

    return  [
        *Xs, 
        c[1],*ve,
        c[2],*vh,
                  ]
pga.graph(graph_func,animate=True,lineWidth=4,grid=True)

```

```python
E,H
```

```python



Xs = map(bm, xs) 
T = lambda x: exp(1/2*(x^P.e0))
Ts = [T(x) for x in Xs]
o = e0.dual()

points = [T>>o for T in Ts]


cga.graph(*points,
          
          conformal=True)

```

```python

pga = Algebra(3,0,1,start_index=0)
P = cga.blades

blade_map_list = [
    (D['e0'], P['e4']),
    (D['e1'], C['e1']),
    (D['e2'], C['e2']),
    (D['e3'], C['e3']),
    (D['e01'], C['e14']),
    (D['e02'], C['e24']),
    (D['e03'], C['e34']),
    (D['e12'], C['e12']),
    (D['e13'], C['e13']),
    (D['e23'], C['e23']),
    (D['e012'], C['e124']),
    (D['e013'], C['e134']),
    (D['e023'], C['e234']),
    (D['e123'], C['e123']),
    (D['e0123'], C['e1234']),
]

bm = BladeMap(blade_map_list)
ni = C.e4+C.e3
no = 0.5*(C.e4-C.e3)
#point = lambda x,y: no + x*C.e1 + y*C.e2 + 0.5*(x*x+y*y)*ni
Xs = map(bm, xs)
point  = lambda x: no + x + 0.5*x**2*ni
T = lambda x: exp(1/2*(x^ni))
Ts = [T(x) for x in Xs]


points = [T>>no for T in Ts]


cga.graph(*points,
          
          conformal=True)

```

```python
BladeMap
```

## Numerical Tests

```python
from kingdon import Algebra
from kingdon.numerical import is_close,exp
from kingdon.calculus import d


D = Algebra(1,3,start_index=0)
locals().update(D.blades)
I  = D.pseudoscalar([1])

K  = D.random_trivector() 
x  = D.random_vector()
F0 = D.random_bivector()

k = K.dual()

assert is_close(I*(k|x) , x^K)
assert is_close(d(D,lambda x: exp(x^K))(x),K*exp(x^K))

assert is_close(exp(K*x) , exp(K|x) *exp(K^x))
assert is_close(exp(K^x), exp((K*x-x*K)/2)) 
assert is_close(exp(K^x), exp(K*x/2)*exp(-x*K/2)) 
assert is_close(exp(x*K)*F0*exp(-K*x) , exp(2*x^K)*exp(x|K)*F0*exp(-x|K))

assert is_close(d(D,lambda x: x^K)(x),K)
assert is_close(d(D,lambda x: x|K)(x),3*K)

```

```python

```

```python
Fx = lambda x: exp(x*K)*F0*exp(-K*x)
Fx(x) *K
```

```python
exp(P^x), exp(P*x/2)*exp(-x*P/2)
```

```python
from kingdon import Algebra 
import numpy as np 
exp = lambda X: X.exp()
normalize = lambda x:(x/np.sqrt((x**2).e))
#cga_norm = lambda x: x/(x|ni).e

def b2pp(B):
    F = normalize(B)
    P  =  0.5*F + 0.5
    P_ = -0.5*F + 0.5
    A = -(P_ * (B | ni))
    B =  (P  * (B | ni))
    return A, B



cga = Algebra(3,1,0,start_index=1)
locals().update(cga.blades)
I = cga.pseudoscalar()
ni = e4+e3
no = 0.5*(e4-e3)
point = lambda x,y: no + x*e1 + y*e2 + 0.5*(x*x+y*y)*ni


T = lambda x1,x2: exp(1/2*ni*(x1*e1 + x2*e2))
R = lambda theta: exp(-1/2*theta*e12)

B1 = no^point(1,0) 
B2 = no^point(0,1) 


cga.graph(b2pp(T(1,0)>>B1),
          grid=True,conformal=True,lineWidth=4)
```

```python
B1*B2
```

viewing field bivectors in  cga , does it help at all.
well, what does this even mean?  the plane reprents all null points, so the idea of a location being a point makes this limited to fields on the null cone 

```python
from kingdon import Algebra 
import numpy as np 
from kingdon.numerical import exp


normalize = lambda x:(x/np.sqrt((x**2).e))
#cga_norm = lambda x: x/(x|ni).e

def bivector_to_point_pair(B):
    F = normalize(B)
    P  =  0.5*F + 0.5
    P_ = -0.5*F + 0.5
    A = -(P_ * (B | ni))
    B =  (P  * (B | ni))
    return A, B

format_for_ganja = lambda x: list(zip( *( [k.flatten() for k in x] ) )) 
 


sta = Algebra(3,1,0,start_index=1)
locals().update(sta.blades)
I = sta.pseudoscalar()
ni = e4+e3
no = 0.5*(e4-e3)
point = lambda x,y: no + x*e1 + y*e2 + 0.5*(x*x+y*y)*ni


dists = np.linspace(1,2,3)
#X1,X2,X3,X4 = np.meshgrid(dists,dists,dists,dists)
#X = X1*e1 + X2*e2 + X3*e3 #+ X4*e4

#P = e123+e124*.33
#F = exp(x^P,30)*(e13)

#[k^P for k in x.flatten()]
X = [sta.vector([x1,x2,x3,0]) for x1 in dists for x2 in dists for x3 in dists]
F = [exp((x|e3)*I)*(e12+.1*e34)  for x in X] 
#As,Bs = bivector_to_point_pair(F)


```

```python
sta.vector([1,2,3,4])
```

```python
sta.graph(*list(zip(As.flatten(), Bs.flatten())),
          grid=True,conformal=True,lineWidth=3)
```

```python
bivector_to_point_pair(F)
```

```python
a=point(1,0)
b=point(1,.5)
B = a^b

normalize = lambda x:(x/np.sqrt((x**2).e))
#cga_norm = lambda x: x/(x|ni).e

def bivector_to_point_pair(B):
    F = normalize(B)
    P  =  0.5*F + 0.5
    P_ = -0.5*F + 0.5
    A = -(P_ * (B | ni))
    B =  (P  * (B | ni))
    return A, B

sta.graph(*bivector_to_point_pair(B),
          bivector_to_point_pair(B),grid=True,conformal=True,lineWidth=3)

```

```python

```

```python

```

```python

```

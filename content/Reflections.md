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

<div style="text-align: right;">
written 02/02/2026<br>
<a href="mailto:alex@810lab.com">alex@810lab.com</a>
</div>

<!-- #region -->
# Reflections

## Summary 
In this we entertain some ideas about reflection coeficient  in special relativity 

## Field decompositions

We are using [Spacetime Algebra](https://en.wikipedia.org/wiki/Spacetime_algebra) , but  with time = $\gamma_4$,  signature = (-1,-1,-1,1)


\begin{align*}
-\gamma_{\not{4}}^2 = \gamma_4^2 = 1 
\end{align*}
Where we use slash to mean  rejection from a given direction. 
<!-- #endregion -->

<!-- #region -->


###  Time
Any object can be it can be decomposed based on symmetries with respect to a given direction. Commonly this is done with respect to time. Given a faraday bivector $F$, and some time direction $\gamma_4$ ,  this leads to a electric/magnetic (E/H) decomp 
\begin{align*}
F &= \frac{1}{2} (F+\gamma_4 F\gamma_4) +\frac{1}{2}( F-\gamma_4F\gamma_4) \\
  &= F\wedge \gamma_4\gamma_4^{-1} +  F\cdot \gamma_4\gamma_4^{-1}\\
  &= F_4 + F_{\not{4}} \\
  &= E+H
\end{align*}

<!-- #endregion -->

<!-- #region -->


 
### General Direction
If instead  of time,  we choose a random direction  $p$ as a decomposition vector, 
\begin{align*}
F &= \frac{1}{2} (F+p Fp) +\frac{1}{2}( F-pFp) \\
  &= F_p + F_{\not{p}} \\
\end{align*}

If we further assume $F$ is given by the plane-wave formula ,  
$$F(x)=e^{p\cdot x I }F_0.$$
where:
* $x$ is a spacetime position,
*  $p=p_4+p_{\not{4}}\equiv k+w$ is the spacetime  propagation vector, and 
* $I \equiv \gamma_{1234}$ is psuedo-scalar. 


then this decomposition is interpreted as a forward and reverse propagating waves. 
$$ F =  F_p + F_{\not{p}} = A+B$$

 (Note, non-standard plane-wave formulas considered [here](PlaneWaves.html))

The decomposition is a choice that is up to us. Here is a picture in $\gamma_{134}$.

![Time Harmonic Normal Incidence](img/FplaneWaveAB_gamma4.svg)


<!-- #endregion -->

## Implementation

```python
from kingdon import Algebra
from kingdon.numerical import exp
sta = Algebra(signature=[-1,-1,-1,1],start_index=1)
locals().update(sta.blades)
F = sta.random_bivector()
E,H  =  1/2*(F + e4*F*e4), 1/2*(F - e4*F*e4)



def decomp(F, p):
    return 1/2*(F + p*F*p), 1/2*(F - p*F*p)

def f(x, p, F0):
    I  = e1234
    return exp((p|x)*I)*F0 

x   = sta.random_vector()
p   = e3 + e4
F0  = e14 + .2*e13  # Ex + Hy
F   = lambda x: f(x,p=p, F0=F0) 
E,H =  decomp(F(x), e4)
A,B =  decomp(F(x), p)
w,k =  decomp(p,e4)
E,H,A,B, E/H, B/A

```

<!-- #region -->


Maxwell's equation in terms of the vector potential $a$
$$\nabla a = F$$
(Sorry for the symbol, we already used $A$ and its a vectory anyway. Also, no one is going to read this anyway )

<!-- #endregion -->

### Simulate it
A EM-field  simulator was made in [this notebook](PlaneWaves.html), below we use it to look  at a simple plane wave with circular polarization. 

```python
#from  fields import * 
#pga.graph(make_fields(F, e12_N=10,e0_bounds=3),  
#          grid=False,lineWidth=4,animate=True,scale=.35,height='300px')
```

<!-- #region -->

## Reflection 
### Time Harmonic Normal Incidence
Electrical Engineering is developed assuming time-harmonic-ness from the outset. Geometrically, this imposes projective geometry in a 'frequency space'.
Consider  a plane wave in a vaccum normally incident upon a surface,  as depicted by the wave-vector diagram below. 


![Time Harmonic Normal Incidence](img/timeHarmonicNormalIncidence_a.svg)

The situtation is time-invariant, and thus $w$ is unchanging. 

The term time-harmonic is an application of [bloch's theorem](https://en.wikipedia.org/wiki/Bloch%27s_theorem#) along time. 

<!-- #endregion -->

## Splits 
### Phase velocity 
$$V_p \equiv wk^{-1} = \frac{w}{k}=\frac{p_4}{p_{\not 4}}= \frac{p \cdot \gamma_4 }{p \wedge\gamma_4 } $$
### Impedance, $Z$
$$Z \equiv EH^{-1}  = \frac{F_4}{F_{\not 4}} = \frac{F\wedge\gamma_4  }{ F\cdot \gamma_4 }$$
### Reflection Coefficient, $\Gamma$
$$\Gamma \equiv AB^{-1}  = \frac{F_p}{F_{\not p}}=  \frac{F\wedge p  }{ F\cdot p }$$




```python
def round(mv, tol=1e-6):
    return mv.filter(lambda c: abs(c) > tol)

def split(F,p):
    num,denom = decomp(F,p)
    return num/denom

Gamma = round(split(F(x),p) )
Z0    = round(split(F(x),e4) )
V     = round(split(p,e4 ))
Gamma, Z0, V
```

```python
 p
```

```python
import numpy as np 

norm = lambda x: x/np.sqrt(abs(((x)**2).e))
na = norm(w) + norm(k)
nb = norm(w) - norm(k)
Na = 1/4*(na*nb)
Nb = 1/4*(nb*na)

p*Nb, p*Na
```

```python
na =  w.normalized() +k.normalized()
nb =  w.normalized() -k.normalized()

import numpy as np 


b = np.tanh(np.pi/4 - np.atanh(((w/k)**2).sqrt().e))*nb
b,na
```

```python
k.normalized()


```

```python


 

```

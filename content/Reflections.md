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

## summary 
In this we entertain some ideas about reflection coeficient  in special relativity 

## Field decompositions

We are using [Spacetime Algebra](https://en.wikipedia.org/wiki/Spacetime_algebra) , but  with time  $\gamma_4$,  signature = (-1,-1,-1,1)

\begin{align*}
\gamma_1^2 = \gamma_2^2 = \gamma_3^2 = -\gamma_4^2 = -1 
\end{align*}

###  Time
Any object can be it can be decomposed based on symmtries with respect to a given direction, which is interpreted as *projection* and *rejection*. Commonly this is done with respect to time. Given a faraday bivector $F$, and some time direction $\gamma_4$ ,  this leads to a electric/magnetic (E/H) decomp 
\begin{align*}
F &= \frac{1}{2} (F+\gamma_4 F\gamma_4) +\frac{1}{2}( F-\gamma_4F\gamma_4) \\
  &= F_4 + F_{\not{4}} \\
  &= E+H
\end{align*}


 
### General Direction
If instead  of time,  we choose a random direction  $p$ as a decomposition vector, 
\begin{align*}
F &= \frac{1}{2} (F+p Fp) +\frac{1}{2}( F-pFp) \\
  &= F_p + F_{\not{p}} \\
\end{align*}

If we further assume $F$ is given by the planewave formula (non-standard plane-wave formulas considered [here](PlaneWaves.html)),  
$$F(x)=e^{p\cdot x I }F_0.$$
where:
* $x$ is a spactime position,
*  $p=p_4+p_{\not{4}}\equiv k+w$ is the spacetime  propagation vector, and 
* $I \equiv \gamma_{1234}$ is psuedo-scalar. 


then this decomposition is interpreted as a forward and reverse propgating waves.
$$ F =  F_p + F_{\not{p}} = A+B$$

 

Here is a picture in $\gamma_{134}$.

![Time Harmonic Normal Incidence](img/FplaneWaveAB_gamma4.svg)

Consert Maxwell's equation in terms of the vector potential $a$
$$\nabla a = F$$
(Sorry for the symbol, we already used $A$ and its a vectory anyway. Also, no one is going to read this anyway )

<!-- #endregion -->

### Simulate it
A EM-field  simulator was made in [this notebook](PlaneWaves.html), below we use it to look  at a simple plane wave with circular polarization. 

```python
from fields import * 

def F(x):
    I  = D.e0123
    F0 = D.e01 + D.e13 # Ex + Hy 
    p  = D.e3 + D.e0   # k + w 
    return exp((p|x)*I)*F0 

pga.graph(make_fields(F, e12_N=10,e0_bounds=3),  grid=False,lineWidth=4,animate=True,scale=.35,height='300px')
```

## Reflection 
### Time Harmonic Normal Incidence
Consider normal incidence  of a plane wave as depicted by the diagram below
![Time Harmonic Normal Incidence](img/timeHarmonicNormalIncidence_a.svg)


## Splits 
### Impedance, $Z$
$$Z \equiv EH^{-1}  = \frac{E}{H}= \frac{F + \gamma_4F\gamma_4 }{ F-\gamma_4F\gamma_4 }$$
### Reflection Coefficient, $\Gamma$
$$\Gamma \equiv AB^{-1}  = \frac{A}{B}= \frac{F + kFk }{ F-kFk }$$





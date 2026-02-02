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
written 01/19/2026<br>
<a href="mailto:alex@810lab.com">alex@810lab.com</a>

</div>



# Maxwell's  Operator

## Summary 
We explore  an operator related to Maxwell's equations, and try to respect differential geometry. 

## Idea 
Consider the following operator :

$$e^{\nabla} = 1+\nabla + \nabla^2 + \ldots$$

Ignore the higher order terms for now. Apply it to the vector potential, $A$.

\begin{align}
e^{\nabla}A &= A+\nabla A + \nabla^2 A = A+F + J
\end{align}

Which is interesting enough to look more at this operator.  Given a definition of the vector derivative in terms of a basis[1], 

$$\nabla \equiv \sum e^i \partial_i .$$
where $e^i \equiv e_i^{-1}$, and $\partial_i \equiv \partial_{e_i}$ is the differential in the $e_i$ direction. This implies 
$$ e^{\nabla} = e^{ \sum e^i \partial_i}  \stackrel{?}{=} \prod e^{ e^i \partial_i}$$

Lets reverse the order of logic and assume $\prod e^{ e^i \partial_i}$ is true. Then under what circumstance does the equality hold:
$$ \prod e^{ e^i \partial_i}  \stackrel{?}{=} e^{ \sum e^i \partial_i}=e^{\nabla}  $$
expanded in a cartesian frame in STA,
$$ e^{  e^0 \partial_0}e^{  e^1 \partial_1}e^{  e^2 \partial_2}e^{  e^3 \partial_3} 
 \stackrel{?}{=} e^{  e^0 \partial_0 + e^1 \partial_1+  e^2 \partial_2+  e^3 \partial_3}  $$

This is almost identical to the transmission line model used in [2]. In the transmission line case, the  equality was justified by the fact that the rotations were small (and in null bivectors).  However, the rotations relative magnitude is what produced the range of behavoir.




## Curvature
For a *linear* vector manifold, the exponential of the differential, $e^{\tau \partial_a}$, implements a translation in the independent variable when operating on a function $f$. 
$$e^{\tau \partial_a}f = f(x+\tau a)    $$

for a small $\tau$. Since translations commute, 

$$ f(x+\tau a +\epsilon b) = 
e^{\tau \partial_a+\epsilon \partial_{b}}\stackrel{?}{=} 
e^{\tau \partial_a}e^{\epsilon \partial_{b}} $$

Again,  reverse the logic and assume the product of exponentials is true in general, 

$$e^{\tau \partial_a}e^{\epsilon \partial_{b}} \stackrel{?}{=} 
e^{\tau \partial_a+\epsilon \partial_{b}} $$

This is true when $\tau$ and $\epsilon$ are small, but small compared to what? Compared to the curvature. This is  what we are assuming by a *linear* vector manifold. One way to measure curvature, would be to move along $a$ then $b$, then back along $-a$ *then* $-b$. Basically reverse your steps in a different order  to get back to the start.  
$$\Delta f 
= e^{-\tau \partial_a}e^{-\epsilon \partial_b}e^{\tau \partial_a}e^{\epsilon \partial_b}f -f
=( e^{-\tau \partial_a}e^{-\epsilon \partial_b}e^{\tau \partial_a}e^{\epsilon \partial_b}-1)f
= (R\tilde{R} -1)f $$
Another method is to compare the difference in $f$ at some nearby point, going either way . 

$$\Delta f 
= e^{\tau \partial_a}e^{\epsilon \partial_b}f -e^{\epsilon \partial_b}e^{\tau \partial_a}f 
= (R -\tilde R )f
$$
Lets see if we can re-work this in terms of conjugation, because conjugation is  dank. 
 $$ e^{\tau \partial_a}e^{\epsilon \partial_b} fe^{\tilde{\epsilon \partial_b}}e^{\tilde{\tau \partial_a}} 
  =f(x+\tau a +\epsilon b  -\epsilon b -\tau a ) =f(x)  $$
which is a dumb thing to write when we are in a linear space. 

<!-- #region -->
### Scalar Curvature 
This section discusses scalar curvature using a non-standard convention for the differential. it is likely not of use. 


The forward differential in the a-direction $ \partial_a $ applied to a function $f$ over a linear manifold is defined by 

\begin{align}
\partial_a f &\equiv \frac{f(x+\tau a )- f(x)}{\tau} \\
\end{align}

Consider defining the backward differential with a reverse symmetry,  
\begin{align}
\  f \partial_a &\equiv \frac{f(x)-f(x-\tau a )}{\tau} \\
\end{align}

The central differential is then seen as an inner product

\begin{align}
\partial_a \cdot f \equiv \frac{1}{2}(\partial_a f + f\partial_a ) 
=\frac{ f(x+\tau a ) - f(x-\tau a )}{2\tau}
\end{align}
The second order central differential is

\begin{align}
\frac{1}{2\tau}(\partial_a f - f\partial_a ) = \frac{1}{\tau}\partial_a \wedge f 

\end{align}

what is 
\begin{align}
\partial_a f \partial_a \\
\end{align}

 the conjugation gives  

\begin{align}
e^{\tau \partial_a} fe^{-\tau \partial_a} &=T_{\tau a}f \tilde T_{\tau a} 

\end{align}

### Exponetial of a  Differential

\begin{align}
\partial_a f &\equiv \frac{f(x+\tau a )- f(x)}{\tau} \\
\partial_a f&= \frac{1}{\tau}(T_{\tau a}-1)f \\
1+ \tau \partial_a  &= T_{\tau a} \\
e^{\tau \partial_a} &=T_{\tau a} 
\end{align}

The backward differential $  \partial_a $  is  
\begin{align}
\  f \partial_a &\equiv \frac{f(x)-f(x-\tau a )}{\tau} \\
f \partial_a &= f (1-T_{-\tau a})\frac{1}{\tau} \\
1- \tau \partial_a  &=T_{-\tau a}\\
e^{-\tau \partial_a} &=\tilde T_{\tau a} 
\end{align}

<!-- #endregion -->

```python

```

## Ortho-not-Normal


## Curvature  on a sphere

```python
from kingdon import Algebra, calculus
from kingdon.numerical import exp  as exp_
exp = lambda x: exp_(x,n=50) # precision control of exp

pga = Algebra(3,0,1,start_index=0) 
locals().update(pga.blades)

f = lambda x:3*x
x= pga.random_vector()
f(x)
lambda f,a: f+da(f,a)

```

## Transmission line Model 

\begin{align}
R &\equiv e_{30} − e_{13} \\
X &\equiv e_{12} − e_{20} \\
G &\equiv e_{30} + e_{13} \\
B &\equiv e_{12} + e_{20}
\end{align}




## Others 
How can split this up by div/curl:

\begin{align}
e^{\nabla}\cdot A &= A + J \\
e^{\nabla}\wedge A &= F
\end{align}




## References

[1]  Hestenes, D. (1984). *Clifford Algebra to Geometric Calculus: A Practical Language for Physics*. D. Reidel Publishing Company.

[2]  A. Arsenovic, "Applications of Conformal Geometric Algebra to Transmission Line Theory," in IEEE Access, vol. 5, pp. 19920-19941, 2017, doi: 10.1109/ACCESS.2017.2727819.


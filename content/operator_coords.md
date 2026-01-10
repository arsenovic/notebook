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
written 01/10/2026<br>
<a href="mailto:alex@810lab.com">alex@810lab.com</a>
</div>


<!-- #region -->
# Operationalizing  Coordinates

Coordinates are modeled as operators acting on the origin in PGA.  

## Context 

When assuming a coordinate system this step is usually glossed over, 
$$ x = g(o)$$

where '$x$' is some position vector (or blade), '$o$' is the origin, and 'g' is some operator.  
A function $f$ of '$x$' is then, $$f(x) = f(g(o))$$
If you do calculus, you have obey the chain rule, and the '$g$'  matters. 

## 2D Cartesian 

The simplest place to start is with a cartesian grid.  to operationalize this, we use [Plane Based Geometric Algebra (PGA)](https://en.wikipedia.org/wiki/Plane-based_geometric_algebra) since it implmenents translations. 



A point is modeled as a dual vector, 

$$\text{origin} = o = e_0^* $$
and a translation by a rotation in a null plane. 
$$T_x = e^{\frac{1}{2}e_0\wedge x}$$
so we have the 
 
$$x = f(o) =  T_x o \tilde{T_x} =  e^{\frac{1}{2}e_0\wedge x}x e^{-\frac{1}{2}e_0\wedge x}$$

If this seems like overkill  then we agree! But adoption this form from the outset allows the other coordinate systems progress more intuitively. We are 'paying it forward'.



Below we implement this with kingdon,  note the operator `>>` implements conjugation; ` T>>o`$= To\tilde{T}$
<!-- #endregion -->

```python
import numpy as np
from kingdon import Algebra
exp = lambda x: x.exp()

pga = Algebra(2,0,1,start_index=0) 
locals().update(pga.blades)

T     = lambda x1,x2: exp(-(x1*e01 + x2*e02)/2)  # translation operator 
dists = np.linspace(-1,1,11)                 # parameters
o     = e0.dual()                            # origin 

points   = [T(x1,x2)>>o for x1 in dists for x2 in dists] # Create points by moving origin around 

pga.graph(*points, grid=False ) 
```

### Add a Grid 
expand the above a little bit by connecting the points together using tuples of points 

```python
import numpy as np
from numpy import pi 
from kingdon import Algebra
exp = lambda x: x.exp()
from itertools import pairwise,chain
c = [0, 1810039, 14245634, 7696563, 15149450, 6727198, 15117058, 10909213, 6710886] 

pga = Algebra(2,0,1,start_index=0) 
locals().update(pga.blades)

T     = lambda x,y: exp(-(x*e01 + y*e02)/2)  # translation operator 
dists = np.linspace(-1,1,11)   # parameters
o     = e0.dual() # origin 
 
 # Create points by moving origin around 
points   = [T(x,y)>>o for x in dists for y in dists]

# create grid by making line-segments connecting points (tuples in ganja)
x_lines = [pairwise( [T(x,y)>>o for y in dists]) for x in dists]
y_lines = [pairwise( [T(x,y)>>o for x in dists]) for y in dists]

pga.graph(
    c[0], *points,
    c[1], *list(chain(*x_lines)),
    c[2], *list(chain(*y_lines)),
    grid=False, lineWidth=3, 
)
```

## 2D Polar 
Next we look at a polar coordinates system which can be implemented with the combination of a Translation and Rotation. 

$$R_\theta = e^{-\frac{\theta}{2} e_{12}}$$

where $\rho$ is the radial distance and $\theta$ is the angle. The combined operator is:

$$RT_{\rho,\theta} = R_\theta T_\rho = e^{-\frac{\theta}{2} e_{12}} e^{-\frac{\rho}{2} e_{01}}$$

This gives us a point in polar coordinates as:

$$x = RT_{\rho,\theta} >> o$$

Below we visualize the polar grid with both angular and radial lines.

```python
pga = Algebra(2,0,1,start_index=0) 
locals().update(pga.blades) 

# operators
T  = lambda x,y: exp(-(x*e01 + y*e02)/2) 
R  = lambda theta: exp(-theta*e12/2)
RT = lambda rho,theta: R(theta)*T(rho,0) 
 
# parameters
angles   = np.linspace(-pi,pi,41)
dists    = np.linspace(.1,2,11)

o = e0.dual() # origin 

 # Create points by moving origin around 
points  = [RT(rho, theta)>>o for theta in angles for rho in dists]
# create grid by connecting points 
angular_lines = [pairwise( [RT(rho, theta)>>o for theta in angles]) for rho in dists]
radial_lines  = [pairwise( [RT(rho, theta)>>o for rho in dists]) for theta in angles]

pga.graph(
    c[0], *points,
    c[1], *list(chain(*angular_lines)),
    c[2], *list(chain(*radial_lines)),
    grid=False, lineWidth=3, 
)
```

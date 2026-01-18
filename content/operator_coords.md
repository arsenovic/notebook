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

## Summary
Coordinates are modeled as operators acting on the origin. This simplifies computing tangent frames, derivatives, etc.

## Context 

When assuming a coordinate system this step is usually glossed over, 
$$ x = g(o)$$

where '$x$' is some position vector (or blade), '$o$' is the origin, and 'g' is some operator.  
A function $f$ of '$x$' is then, $$f(x) = f(g(o))$$
If you do calculus, you have obey the chain rule, and the '$g$'  matters. 

## 2D Cartesian 

The simplest place to start is with a 2D cartesian coordinate system.  to operationalize this, we use [Plane Based Geometric Algebra (PGA)](https://en.wikipedia.org/wiki/Plane-based_geometric_algebra) since it implmenents translations. 



A point is modeled as a dual vector, 

$$\text{origin} = o = e_0^* $$
and a translation by a rotation in a null plane. 
$$T_x = e^{\frac{1}{2}e_0\wedge x}$$

so we have our position function, f 
 
$$x = f(o) =  T_x o \tilde{T_x} =  e^{\frac{1}{2}e_0\wedge x}x e^{-\frac{1}{2}e_0\wedge x}$$

If this seems like overkill  then we agree! But adopting an operator perspective from the outset allows the other coordinate systems progress more intuitively. You are 'paying it forward'.

We implement a function that operates over scalar parameters (aka coordinates) $x_1, x_2$, so the literal opertor will be 

$$T_{x_1, x_2}=  e^{\frac{1}{2}e_0\wedge(x_1 e_1 + x_2e_2)}$$



Below we implement this with kingdon,  note the operator `>>` implements conjugation; ` T>>o`$= To\tilde{T}$
<!-- #endregion -->

```python
import numpy as np
from kingdon import Algebra
exp = lambda x: x.exp()

pga = Algebra(2,0,1,start_index=0) 
locals().update(pga.blades)

T     = lambda x1,x2: exp(-e0*(x1*e1 + x2*e2)/2)  # translation operator 
dists = np.linspace(-1,1,11)                 # parameters
o     = e0.dual()                            # origin 

points   = [T(x1,x2)>>o for x1 in dists for x2 in dists] # Create points by moving origin around 
  

pga.graph(*points, grid=False ) 
```

### Add a Grid 
expand the above a little bit by connecting the points together using tuples of points 

```python
from itertools import pairwise,chain
 
# create grid by making line-segments connecting points (tuples in ganja)
h_lines = [pairwise( [T(x1,x2)>>o for x1 in dists]) for x2 in dists]
v_lines = [pairwise( [T(x1,x2)>>o for x2 in dists]) for x1 in dists]


pga.graph(
    *list(chain(*h_lines)),
    *list(chain(*v_lines)),
    *points,
    grid=False, lineWidth=3, 
)
```

## 2D Polar 
Next we look at a polar coordinates system which can be implemented with the combination of a Translation and Rotation. 

$$R_\rho = e^{-\frac{\rho}{2} e_{01}}, \qquad R_\theta = e^{-\frac{\theta}{2} e_{12}}$$

where $\rho$ is the radial distance and $\theta$ is the angle. The combined operator is:

$$RT_{\rho,\theta} = R_\theta T_\rho = e^{-\frac{\theta}{2} e_{12}} e^{-\frac{\rho}{2} e_{01}}$$

This gives us a point in polar coordinates as:

$$x =  e^{-\frac{\theta}{2} e_{12}} e^{-\frac{\rho}{2} e_{01}}o   e^{\frac{\rho}{2} e_{01}}e^{\frac{\theta}{2} e_{12}} ,$$

which we implement with  $$ RT_{\rho,\theta} >> o$$

Below we visualize the polar grid with both angular and radial lines.

```python
from numpy import pi 

# operators
T  = lambda x1,x2: exp(-e0*(x1*e1 + x2*e2)/2) 
R  = lambda theta: exp(-theta*e12/2)
RT = lambda rho,theta: R(theta)*T(rho,0) 
 
# parameters
angles   = np.linspace(-pi,  pi,41)
dists    = np.linspace(.1,2,11)

o = e0.dual() # origin 

 # Create points by moving origin around 
points  = [RT(rho, theta)>>o for theta in angles for rho in dists]
# create grid by connecting points 
angular_lines = [pairwise( [RT(rho, theta)>>o for theta in angles]) for rho in dists]
radial_lines  = [pairwise( [RT(rho, theta)>>o for rho in dists]) for theta in angles]

c = [0, 1810039, 14245634, 7696563, 15149450, 6727198, 15117058, 10909213, 6710886]  # colors 
pga.graph(
    c[1], *list(chain(*angular_lines)),
    c[2], *list(chain(*radial_lines)),
    c[0],  *points,
    grid=False, lineWidth=4, 
)
```

### Symbolics 
kingdon can make use of sympy symbolics. 

```python
from sympy import symbols, sin, cos, pretty
# Use Greek letter unicode in symbol names for display
rho, theta = symbols('ρ θ', real=True, positive=True)

R_  = lambda t: cos(t/2)*e - sin(t/2)*e12
T_  = lambda r: 1 + e01*r/2
RT_ = lambda r, t: R_(t)*T_(r)
o   = e0.dual()
(RT_(rho, theta) >> o).dual().proj(e12)

```

### Tangent Frame 
Since we have operationalized the coordinate function, computing a tangent frame is simple. 

```python
# create a tangent frame at each point (vectors are tuples of 2 points)
tmag = .1 # tangent vector magnitude
radial   = (o, T(tmag,0)>>o) # tangent vector in radial direction
angular  = (o, T(0,tmag)>>o) # tangent vector in angular direction
radials  = [RT(rho, theta)>>radial for theta in angles for rho in dists]
angulars = [RT(rho, theta)>>angular for theta in angles for rho in dists]

pga.graph(
    c[1], *radials,
    c[2], *angulars,
    c[0],  *points,
    grid=False, lineWidth=4, )

```

## 3D Polar (Latitude, Longitude)

Extending to 3D spherical coordinates, we combine two rotations with a radial translation. The latitude $\phi$ and longitude $\theta$ rotations are:

$$R_\theta = e^{-\frac{\theta}{2} e_{12}}, \quad R_\phi = e^{-\frac{\phi}{2} e_{13}}$$

The combined operator for a point at radius $\rho$ is:

$$RT_{\rho,\theta,\phi} = R_\phi R_\theta T_\rho = e^{-\frac{\phi}{2} e_{13}} e^{-\frac{\theta}{2} e_{12}} e^{-\frac{\rho}{2} e_{01}}$$

"Null Island" is the random reference point  on earth where latitude and longitude are both zero.  

Below we visualize the globe with latitude and longitude grid lines.

```python
pga = Algebra(3,0,1,start_index=0) 
locals().update(pga.blades)

# operators
T  = lambda x1,x2,x3: exp(-e0*(x1*e1 + x2*e2 + x3*e3)/2) 
R  = lambda theta, phi : exp(-theta*e12/2) * exp(-phi*e13/2) 
RT = lambda rho,theta,phi: R(theta,phi)*T(rho,0,0) 

# parameters
thetas = np.linspace(-pi,pi,21)  
phis   = np.linspace(-pi/2,pi/2,11)[1:-1]  
rho    = 1 # earth radius

o           = e0.dual() # origin 
null_island = T(rho,0,0)>>o # lat=lon=0 point 
points      = [R(theta, phi)>>null_island for theta in thetas for phi in phis ]

# create grid by connecting points 
lat_lines = [pairwise( [R(theta,phi)>>null_island for theta in thetas]) for phi in phis]
lon_lines = [pairwise( [R(theta,phi)>>null_island for phi in phis]) for theta in thetas]

pga.graph(
    c[0], *points,
    c[5], null_island,
    c[1], *list(chain(*lat_lines)),
    c[2], *list(chain(*lon_lines)),
    grid=False, lineWidth=3,
)
```

### Tangent Frame 
next to compute the tangent frame. In geodesy, the standard tangenth frame is defined as locally [East, North, Up]'.

We also vectorize the parameters using meshgrid.

```python
# def operators 
T  = lambda x1,x2,x3: exp(-e0*(x1*e1 + x2*e2 + x3*e3)/2) 
R  = lambda theta, phi : exp(-theta*e12/2)*exp(-phi*e13/2) 

# parameters
thetas = np.linspace(-pi,pi,21)              
phis   = np.linspace(-pi/2,pi/2,11)[1:-1]  
#theta, phi = np.meshgrid(thetas,phis)
rho    = 1.5

# objects
o     = e0.dual() # origin 
tmag  = rho/10 # magnitude of tangent vectors
east  = (o, T(0,tmag,0)>>o) # e2
north = (o, T(0,0,tmag)>>o) # e3
up    = (o, T(tmag,0,0)>>o) # e1

# create operator
Rs     = [R(theta,phi)*T(rho,0,0) for theta in thetas for phi in phis]

# operate 
points = [R>>o for R in Rs]
easts  = [R>>east for R in Rs]
norths = [R>>north for R in Rs]
ups    = [R>>up for R in Rs]
 
pga.graph(
    c[0], *points,
    c[1], *easts,
    c[2], *norths,
    c[3], *ups,
    grid=False, lineWidth=5,
)
```

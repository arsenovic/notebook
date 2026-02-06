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
written 01/17/2026<br>
<a href="mailto:alex@810lab.com">alex@810lab.com</a>
</div>


# Plane Waves 

## Summary 
Starting from the electromagnetic plane wave equation in [Spacetime Algebra](https://en.wikipedia.org/wiki/Spacetime_algebra), a similar equation is reverse engineered which follows conjugation. It is not known to be useful. 

<!-- #region -->
## Idea

Consider  the form of a EM plane wave given in Hestenes' Spacetime Algebra 

$$  F_x\equiv F(x) = e^{I(k\cdot x)}F_0 $$

where;
* $x$ is a location in spacetime
* $k$ is the propagation vector
* $I$ is psuedoscalar. 

This equation is interpreted as a duality rotation of a fixed bivector $F_0$, parameterized by $k$, with $x$ being  an independent variable.
The exponent can be re-written as 

$$  I(k\cdot x)= (Ik)\wedge x =  K\wedge x$$


where $K$ is a trivector dual to $k$. ( We prefer trivectors  for propagation constants since the units match the grade ( rad/m rad/s), but keeping with dual vectors works just as good.) The  function

$$ F_x=  e^{K\wedge x}F_0$$
 solves the equation, 
$$ \nabla F_x= KF_x. $$

<!-- #endregion -->

So far we have changed nothing from STA. If we re-express $F_x$ in terms of the geometric product, it suggests to us a formula which employs conjugation, such as  

$$ e^{K\wedge x}F_0 = e^{Kx-xK}F_0 \stackrel{?}{=} e^{Kx}F_0e^{-xK}. $$

(we ignore all factors of 1/2 untill we need to precise answer).  For the last  equality to hold, $F_0$ must commute with the bivector $K\cdot x$. However, since $x$ is an independent variable, this cannot be guaranteed without more assumptions. So,  there are two possibilities 
1) This reverse engineered formula is useless
2) The plane-wave forumula is making implicit assumptions which justify the commutation assumption. 

For now, lets proceed and  see what such a forumla provides. Expanding the product 

$$xK = x\cdot K + x \wedge K  = B+\alpha I  $$

Where $B$ is a bivector. Since the psuedoscalar commutes with $F_0$, we can then write, 


$$  e^{ K\wedge x} e^{K\cdot x} F_0e^{-K \cdot x} = e^{ \alpha I }e^{B}F_0e^{-B}   $$

This is interepreted as a STA **position-dependent lorentz rotation** (and duality rotation).  Along certain axis, namely $x\cdot K = 0$, the lorentz transforms disappears and this reduces to the normal plane wave. It is reasonable to extend this to spacially varying $K$, $K\rightarrow K(x)$.



## Simulate
Since is somewhat easier to simulate than to work out the maths, below are some field simulations of this wave function. The simulations generate the field for a variety of $K$'s . The fields are mapped from STA into 3dPGA for visualization. Since kingdown allows pytorch.tensors as coefficient this computation is all done on the GPU. 

<div style="display: flex; justify-content: center; gap: 80px; flex-wrap: wrap;">
  <div>
    <video width="400" autoplay muted loop>
      <source src="https://www.dropbox.com/scl/fi/aomtyu22x4c4f99vx0d0n/traveling_wave.mp4?rlkey=kr37kpdy4kgejnvgcit53au19&st=wnvwnpgy&dl=1&raw=1" type="video/mp4">
      Your browser does not support the video tag.
    </video>
  </div>    

  <div>
    <video width="400" autoplay muted loop>
      <source src="https://www.dropbox.com/scl/fi/i2dr9ka3t5j7i4wh97gdq/oscillator.mp4?rlkey=dgs6rdcrnn96etz5qe77zskhb&st=z2gksssc&dl=1&raw=1" type="video/mp4">
      
      Your browser does not support the video tag.
    </video>
  </div>
</div>

<div style="display: flex; justify-content: center; gap: 80px; flex-wrap: wrap;">
  <div>
    <video width="400" autoplay muted loop>
      <source src="https://www.dropbox.com/scl/fi/5rrkecj9p9bk4ki47nlfv/sheets_wave.mp4?rlkey=qhi28km6v8pwiphrijfrp1fn0&st=yggsxk8l&dl=1&raw=1" type="video/mp4">
      Your browser does not support the video tag.
    </video>
  </div>    

<div style= "">
  <img src="https://www.dropbox.com/scl/fi/ujvinrv8dttc2qfath3wh/P-d123-p5d012-A0-d1-sweep_v3mag.gif?rlkey=9z06knp3j0ddfzou2ccgprr60&st=mgp14a0g&dl=1&raw=1" width="400" alt="Swirl animation">
</div>
</div>

 



## Extensions
### A-field 
The F-field can also be generated from a vector potential field, $A$. The same equation is used for $A$, and then the $F$ is computed by using a numerical implementation of the geometric derivative. 

* trivector pregenerator becomes trivector+vector to model loss $K\rightarrow K+k  \stackrel{?}{=} e^{\alpha I }K$
* extend  $F$ to be a arbitrary multivector $M$
* trivector/vector pre-generator to also be a function of position, $K\rightarrow K(x)$ 

$$ M_x = e^{K_x x}M_0e^{-xK_x}. $$

 ## Details 
 * order of $xK$ vs $Kx$,  and + vs -
 * re-work in terms fo vector $k$ 




## Numerical Tests
Tests of the maths to ensure we are consistent. 

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

# try to run this block up to 3 times, and break if all assertions pass
# need this because its seeded with random values, and its numerical 
for _ in range(3):
    try:
        assert is_close(I*(k|x) , x^K)
        assert is_close(d(D,lambda x: exp(x^K))(x),K*exp(x^K))

        assert is_close(exp(x*K) , exp(x|K) *exp(x^K))
        assert is_close(exp(x^K), exp((x*K-K*x)/2)) 
        assert is_close(exp(x^K), exp(x*K/2)*exp(-K*x/2)) 
        assert is_close(exp(x*K)*F0*exp(-K*x) , exp(2*x^K)*exp(x|K)*F0*exp(-x|K))
        assert is_close(exp(x*K)*F0*exp(-K*x) , exp(x*K)>>F0)

        assert is_close(d(D,lambda x: x^K)(x),K)
        assert is_close(d(D,lambda x: x|K)(x),3*K)
    except AssertionError:
        continue
    break
```

## Field Visualization 
code used to generate the Field-animations. Various values for K, x .

```python
from kingdon import Algebra
from kingdon.numerical import exp  as exp_
exp = lambda x: exp_(x,n=50) # precision control of exp

from kingdon.calculus import d
from kingdon.blademap import BladeMap
import numpy as np 
import torch
from timeit import default_timer
c = [0, 1810039, 14245634, 7696563, 15149450, 6727198, 15117058, 10909213, 6710886]  # colors 
 
  
# Use GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)
## Set up algebras and blade map
sta = Algebra(1,3,start_index=0)
D   = sta.blades
DI  = sta.pseudoscalar([1])

pga = Algebra(3,0,1)
P   = pga.blades
up  = lambda x: (x +  P.e0).broadcast().dual()
to_numpy = lambda x:x.cpu().numpy()

bm  = BladeMap(alg1=sta, alg2 =pga) 
torch_kw = {'device':device, 'dtype':torch.float64} # you need float64 for calculus to work.

def make_x(e12_N, e12_bounds, e0_N, e0_bounds, e3_N, e3_bounds):
    ''' Create a grid of 'x' (events)  '''
    x0_range = torch.linspace(-e0_bounds, e0_bounds, e0_N, **torch_kw)
    x1_range = torch.linspace(-e12_bounds, e12_bounds, e12_N, **torch_kw)
    x2_range = torch.linspace(-e12_bounds, e12_bounds, e12_N, **torch_kw)
    x3_range = torch.linspace(-e3_bounds, e3_bounds, e3_N, **torch_kw)
    x0,x1,x2,x3= torch.meshgrid(x0_range, x1_range, x2_range, x3_range, indexing='ij')
    return x0*D.e0 + x1*D.e1 + x2*D.e2 + x3*D.e3

def render_frames(x,Fx,vector_scale=1):
    '''Create E/H fields representation in PGA from a F-field in STA and
    render each frame for a ganja.js visualization
    Args:
        x: grid of events (vectors)
        Fx: EM field   at each event (bivectors)
        vector_scale: scaling factor for the field vectors'''
    E = Fx.proj(D.e0)
    H = Fx.proj(D.e0.dual).dual()
    frames = []
    N_frames = x.shape[1]

    # TODO move x_, E_, H_ calculation outside loop
    # TODO pass numpy array straight to pga.graph to avoid flattening
    for frame_idx in range(N_frames):
        mid = int(N_frames/2)
        s = (frame_idx, slice(None), slice(None), slice(None))
        x_, E_, H_ = [X[s] for X in (x, E, H)]
        x_ = up(bm(x_.proj(D.e123)))
        E_ = bm(E_)
        H_ = bm(H_)
        
        vmag = 1/x.shape[2]*vector_scale
        x_flat = x_.map(to_numpy).flatten()
        ve_flat = (exp(vmag*E_).broadcast()>>x_).map(to_numpy).flatten()
        vh_flat = (exp(vmag*H_).broadcast()>>x_).map(to_numpy).flatten()
        
        
        frames.append([c[1], *[(a,b) for a,b in zip(x_flat, ve_flat)],
                       c[2], *[(a,b) for a,b in zip(x_flat, vh_flat)]])
    return frames

def make_graph_func(frames,speed ):
    '''light wrapper to create a graph function that cycles through frames'''
    def graph_func():
        t = int(default_timer()*speed)
        num_frames = len(frames)
        
        cycle_length = 2 * (num_frames - 1)
        t_mod = t % cycle_length
        
        if t_mod < num_frames:
            frame_idx = t_mod
        else:
            frame_idx = cycle_length - t_mod
        
        return frames[frame_idx]
    return graph_func


scale=1/(2*2)
speed = 50
vector_scale =2
x = make_x(e12_bounds = 6,
           e0_bounds  = 4,
           e3_bounds  = 0,
           e12_N      = 21,
           e0_N       = 100,
           e3_N       = 1,
           )
def F(x):
    F0 = D.e01 + D.e13
    K  = (D.e123 + .01 * D.e012)#*exp(-.1*D.e0123)
    xK = x*K
    return exp(xK) >> F0

def A(x):
    A0 = D.e3
    #K  = lambda x: exp((.1*(D.e0 )|x)*D.e12)>>(D.e123 +  .5*D.e013)#*exp(-1*D.e0123)
    K  = lambda x: (D.e12/D.e3 +  .5*D.e13/D.e0)#*exp(-1*D.e0123)
    xK = x*K(x)
    return exp(xK)>>A0

# Fx = F(x) # make F-field directly
Fx =  d(sta, A)(x)  # make F-field from potential A
frames     = render_frames(x, Fx, vector_scale=vector_scale)
graph_func = make_graph_func(frames, speed)
vector_scale = float(Fx[0].values()[0].max())

## uncomment to animate
pga.graph(graph_func, grid=True,lineWidth=2,animate=True,scale=.15,height='700px')

```

And thats all we got for now. 

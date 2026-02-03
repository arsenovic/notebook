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

```python
from kingdon import Algebra
from kingdon.numerical import exp  as exp_
exp = lambda x: exp_(x,n=50) # precision control of exp

from kingdon.calculus import d
from kingdon.blademap import BladeMap
import numpy as np 
import torch
from timeit import default_timer

def blades_in(x,tol=1e-9):
    """
    Returns the list of blade names present in this multivector.
    
    :return: List of blade names as strings, e.g. ['e1', 'e23', 'e123']
    """
    result = []
    for i, k in enumerate(x._keys):
        val = x._values[i]
        # Handle both scalar and array coefficients
        if np.any(np.abs(val) > tol):
            result.append(x.algebra.bin2canon[k])
    return result

c = [0, 1810039, 14245634, 7696563, 15149450, 6727198, 15117058, 10909213, 6710886]  # colors 


device = torch.device("cuda" if torch.cuda.is_available() else "cpu");print("Using device:", device)
torch_kw = {'device':device, 'dtype':torch.float64} # you need float64 for calculus to work.


# e0 is dual origin, e4 = time , e123 is space
stap = Algebra(signature = [0,1,1,1,-1],start_index=0)
pga  = Algebra(3,0,1)
bm   = BladeMap(alg1=stap, alg2 =pga) 

S    = stap.blades
locals().update(S)
I    = stap.pseudoscalar([1])
up   = lambda x: ((x+e0).dual()|e4).broadcast() # pga up() done in stap
to_numpy = lambda x:x.cpu().numpy()


 

et,ek = e4,e3 

def make_x(e12_bounds, e12_N, e3_bounds, e3_N, e4_bounds, e4_N, use_torch=False):
    if use_torch:
        x1_range = torch.linspace(-e12_bounds, e12_bounds, e12_N, **torch_kw)
        x2_range = torch.linspace(-e12_bounds, e12_bounds, e12_N, **torch_kw)
        x3_range = torch.linspace(-e3_bounds, e3_bounds, e3_N, **torch_kw)
        x4_range = torch.linspace(-e4_bounds, e4_bounds, e4_N, **torch_kw)
        x4,x1,x2,x3= torch.meshgrid(x4_range, x1_range, x2_range, x3_range, indexing='ij')
    else:
        x1_range = np.linspace(-e12_bounds, e12_bounds, e12_N)
        x2_range = np.linspace(-e12_bounds, e12_bounds, e12_N)
        x3_range = np.linspace(-e3_bounds, e3_bounds, e3_N)
        x4_range = np.linspace(-e4_bounds, e4_bounds, e4_N)
        x4,x1,x2,x3= np.meshgrid(x4_range, x1_range, x2_range, x3_range, indexing='ij')
    return x4*e4 + x1*e1 + x2*e2 + x3*e3  

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


def F(x):
    F0 = e1*et + e1*ek
    k  = ek+.01*et
    xk = x*k
    return exp(xk) >> F0

x = make_x( # make a indepenent variable 'x'
    e12_bounds=1, e12_N=1, 
    e3_bounds =1,  e3_N=1,
    e4_bounds =1,  e4_N=1, 
    use_torch=False)
 

Fx = F(x)
mag =1 

list(map(blades_in, [E,H,Be,Bh,ve,vh,up(x)]))
blades_in((x+e0).dual()|e4)
```

```python
def render_frames(x,Fx, mag=1):
    E  = (Fx + (e4>>Fx))/2 # E/H field decomposition
    H  = (Fx - (e4>>Fx))/2
    Be = E*e04      # generator for E/H field vectors (ex: e12->e3->(e3+e0).dual(in pga))
    Bh = H*e0123            
    X  = up(x)
    ve = exp(mag*Be)>>X
    vh = exp(mag*Bh)>>X
    
    
    
    frames = []
    N_frames = x.shape[1]
    for frame_idx in range(N_frames):
        s = (frame_idx, slice(None), slice(None), slice(None))
        X_, E_, H_ = [q[s] for q in (X, E, H)]
        x_flat,ve_flat,vh_flat= [bm(k).flatten() for k in (X_, E_, H_)]  # flatten arrays
        frames.append(
            [c[1], *[(a,b) for a,b in zip(x_flat, ve_flat)],
            c[2], *[(a,b) for a,b in zip(x_flat, vh_flat)]])
    return frames

def F(x):
    F0 = e1*et + e1*ek
    k  = ek+.01*et
    xk = x*k
    return exp(xk) >> F0

x = make_x(
    e12_bounds=1, e12_N=20, 
    e3_bounds=1, e3_N=3,
    e4_bounds=1, e4_N=4, 
    use_torch=False)

speed=1

Fx = F(x)
X  = up(x)
frames = render_frames(x,Fx,mag=.10)
graph_func = make_graph_func(frames, speed)
pga.graph(*frames[0], grid=False,lineWidth=2,scale=1,height='700px')

```

```python
x.shape
```

```python

## uncomment to animate
#stap.graph(graph_func, grid=True,lineWidth=2,animate=True,scale=.15,height='700px')
pga.graph(*frames[0], grid=False,lineWidth=2,scale=.15,height='700px')
#pga.graph(*x_flat, grid=False,lineWidth=2,scale=.15,height='700px')

```

```python
frames[0]
```

```python
up(x)
```

```python

```

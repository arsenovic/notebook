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

# Animation parameters
num_frames = 120  # Total number of precomputed frames
animation_duration = 8.0  # Duration of one complete cycle in seconds

# Precompute all frames
N=11
bounds = 4
x1, x2 = np.mgrid[-bounds:bounds:N*1j, -bounds:bounds:N*1j]

precomputed_frames = []

for frame_idx in range(num_frames):
    t = frame_idx / num_frames * 4 * np.pi  # sweep through one full cycle
    
    x = x1*D.e1 + x2*D.e2 + 3*np.sin(t)*D.e0
    xs = x.flatten()

    F0 = D.e01 #+ D.e02.dual # E1+H2
    K = D.e123+.5*D.e012
    #F = lambda x: exp(K*x)>>F0

    #Fs = [F(x) for x in xs]
    #E,H = zip(*[(F.proj(D.e0), F.proj(D.e0.dual).dual()) for F in Fs])
    F= exp(K*x)>>F0
    E = (F.proj(D.e0)).flatten()
    H = (F.proj(D.e0.dual).dual()).flatten()
    
    E = list(map(bm,E))
    H = list(map(bm,H))
    xs_proj = [k.proj(D.e12) for k in xs]
    Xs = list(map(up, (map(bm, xs_proj))))
    vmag = .5*1/N
    ve = [(x,.5*exp(vmag*k)>>x) for x,k in zip(Xs,E)]
    vh = [(x,.5*exp(vmag*k)>>x) for x,k in zip(Xs,H)]

    precomputed_frames.append([ c[1],*ve, c[2],*vh])

print(f"Precomputed {len(precomputed_frames)} frames")
print(f"Animation duration: {animation_duration} seconds per cycle")

def graph_func():
    frame_idx = int((default_timer() / animation_duration * num_frames) % num_frames)
    return precomputed_frames[frame_idx]

pga.graph(graph_func,animate=True,lineWidth=4,grid=False)

```

```python
import cProfile
import pstats
from io import StringIO

# Profile the frame computation
pr = cProfile.Profile()
pr.enable()

# Recompute frames with profiling
num_frames_profile = 10  # Profile on fewer frames for speed
precomputed_frames_profile = []

for frame_idx in range(num_frames_profile):
    t = frame_idx / num_frames_profile * 4 * np.pi
    
    x = x1*D.e1 + x2*D.e2 + 5*np.sin(t)*D.e0
    xs = x.flatten()

    F0 = D.e01 + D.e02.dual
    K = D.e123+.5*D.e012
    F = lambda x: exp(K*x)>>F0

    Fs = [F(x) for x in xs]
    E,H = zip(*[(F.proj(D.e0), F.proj(D.e0.dual).dual()) for F in Fs])
    
    E = list(map(bm,E))
    H = list(map(bm,H))
    xs_proj = [k.proj(D.e12) for k in xs]
    Xs = list(map(up, (map(bm, xs_proj))))
    vmag = .5*1/N
    ve = [(x,.5*exp(vmag*k)>>x) for x,k in zip(Xs,E)]
    vh = [(x,.5*exp(vmag*k)>>x) for x,k in zip(Xs,H)]

pr.disable()

# Print profiling results
s = StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
ps.print_stats(20)  # Top 20 functions
print(s.getvalue())

```

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
    t = default_timer() /2
    #print(t)

    N=15
    bounds = 2
    x1, x2 = np.mgrid[-bounds:bounds:N*1j, -bounds:bounds:N*1j]
    x = x1*D.e1 + x2*D.e2 + 5*np.sin(t)*D.e0
    xs= x.flatten()

    F0 = D.e01 #+ D.e02.dual # E1+H2
    K = D.e123+.1*D.e012
    F = lambda x: exp(K*x)>>F0

    
    F= exp(K*x)>>F0
    E = (F.proj(D.e0)).flatten()
    H = (F.proj(D.e0.dual).dual()).flatten()
    ## convert to pga to display 
    E = list(map(bm,E))
    H = list(map(bm,H))
    xs = [k.proj(D.e12) for k in xs]
    Xs = list(map(up, (map(bm, xs))))
    vmag = .5*1/N
    ve = [(x,.5*exp(vmag*k)>>x) for x,k in zip(Xs,E)]
    vh = [(x,.5*exp(vmag*k)>>x) for x,k in zip(Xs,H)]

    return  [
        #*Xs, 
        c[1],*ve,
        c[2],*vh,
                  ]
pga.graph(graph_func,animate=True,lineWidth=4,grid=False)

```

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


 
t=1
N=7
bounds = 2
x1, x2 = np.mgrid[-bounds:bounds:N*1j, -bounds:bounds:N*1j]
x = x1*D.e1 + x2*D.e2 + 5*np.sin(t)*D.e0
xs= x.flatten()

F0 = D.e01 + D.e02.dual # E1+H2
K = D.e123+.5*D.e012
F= exp(K*x)>>F0
E = (F.proj(D.e0)).flatten()
H = (F.proj(D.e0.dual).dual()).flatten()


E = list(map(bm,E))
H = list(map(bm,H))
xs = [k.proj(D.e12) for k in xs]
Xs = list(map(up, (map(bm, xs))))
vmag = .5*1/N
ve = [(x,.5*exp(vmag*k)>>x) for x,k in zip(Xs,E)]
vh = [(x,.5*exp(vmag*k)>>x) for x,k in zip(Xs,H)]

 
```

```python
F= exp(K*x)>>F0
E = (F.proj(D.e0)).flatten()
H = (F.proj(D.e0.dual).dual()).flatten()

```

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
    t = default_timer() 
    #print(t)

    N=2
    x1, x2 = np.mgrid[-1:1:N*1j, -1:1:N*1j]
    x = x1*D.e1+np.sin(t)*x2*D.e2 
    xs= x.flatten()
 
    Xs = list(map(up, (map(bm, xs))))
 
    return  [
        *Xs, 
 
                  ]
pga.graph(graph_func,animate=True,lineWidth=4,grid=True)

```

```python

```

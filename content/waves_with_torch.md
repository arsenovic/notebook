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
exp = lambda x: exp_(x,n=30)

from kingdon.calculus import d
from kingdon.blademap import BladeMap
import numpy as np 
import torch
from timeit import default_timer
c = [0, 1810039, 14245634, 7696563, 15149450, 6727198, 15117058, 10909213, 6710886]  # colors 



def to_numpy(X):
    """
    Convert a multivector with PyTorch tensor coefficients to NumPy array coefficients.
    
    If coefficients are already NumPy arrays or scalars, they are returned unchanged.
    PyTorch tensors are converted to NumPy via .detach().cpu().numpy().
    
    :return: A new MultiVector with NumPy array coefficients.
    """
    def convert_value(v):
        # Check if it's a PyTorch tensor
        if hasattr(v, 'detach'):
            # PyTorch tensor - convert to NumPy
            return v.detach().cpu().numpy()
        elif isinstance(v, np.ndarray):
            # Already NumPy
            return v
        else:
            # Scalar or other type - leave as is
            return v
    
    converted_values = [convert_value(v) for v in X._values]
    return X.fromkeysvalues(X.algebra, keys=X._keys, values=converted_values)

# Use GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")


sta = Algebra(1,3,start_index=0)
D = sta.blades
DI  = sta.pseudoscalar([1])

N = 30
bounds = 5
speed = 20
vector_scale =2

x0_range = torch.linspace(-3, 3, 60, device=device)
x1_range = torch.linspace(-bounds, bounds, N, device=device)
x2_range = torch.linspace(-bounds, bounds, N, device=device)
x3_range = torch.linspace(-0, bounds, 1, device=device)

x0,x1,x2,x3= torch.meshgrid(x0_range, x1_range, x2_range, x3_range, indexing='ij')
x = x0*D.e0 + x1*D.e1 + x2*D.e2 + x3*D.e3 

# compute teh field  and split into E/Hs
F0 = D.e01+D.e13
K  = D.e123# + 0.1 * D.e012
F = exp(K * x) >> F0
E = F.proj(D.e0)
H = F.proj(D.e0.dual).dual()


pga = Algebra(3,0,1)
P = pga.blades
bm = BladeMap(alg1=sta, alg2 =pga) 
up = lambda x: (x +  P.e0).broadcast().dual()


## convert to pga to display  
def graph_func():
    t = int(default_timer()*speed) % x.shape[1]
    e0_idx = t
    
    e3_idx = 0
    s = (e0_idx, slice(None), slice(None), e3_idx)

    x_, E_,H_ = [X[s] for X in (x,E,H)]


    x_ = up(bm(x_.proj(D.e123))) # throw away time, since we need e0 for projection
    E_ = bm(E_)
    H_ = bm(H_)


    vmag = 1/N*vector_scale

    x_flat =  to_numpy(x_).flatten()
    ve_flat = to_numpy (exp(vmag*E_).broadcast()>>x_).flatten()
    vh_flat = to_numpy (exp(vmag*H_).broadcast()>>x_).flatten()

    return [#*x_flat,
            c[1],*[(a,b) for   a,b in zip(x_flat, ve_flat)],
            c[2],*[(a,b) for   a,b in zip(x_flat, vh_flat)],]
pga.graph(graph_func, grid=False,lineWidth=4,animate=True,scale=.1)
```

```python

```

```python

```

```python

F0 = D.e01
K  = D.e123 + 0.1 * D.e012
 
F = exp(K * x) >> F0
E = F.proj(D.e0)
H = F.proj(D.e0.dual).dual()

pga = Algebra(3,0,1)
P = pga.blades
bm = BladeMap(alg1=sta, alg2 =pga) 
up = lambda x: (x +  P.e0).broadcast().dual()

## convert to pga to display  
sl = lambda X,e0_idx, e3_idx: X[e0_idx][:][:][e3_idx]  #slice out time and one space
e0_idx,e3_idx = 0,-1
x, E,H = [sl(X,e0_idx, e3_idx) for X in (x,E,H)]

x = up(bm(x.proj(D.e123))) # throw away time, since we need e0 for projection
E = bm(E)
H = bm(H)


vmag = .5*1/N

x_flat = x.flatten()
#ve = [(X, Ve) for X,Ve  in zip(x_flat, (exp(vmag*E).broadcast()>>x).flatten())]
pga.graph(*to_numpy(x).flatten(),grid=True)
```

```python
def tensors_to_multivectors(tensor_list, algebra):
    """Convert a list of tensors to a list of MultiVectors in the given algebra."""
    from kingdon.multivector import MultiVector
    result = []
    for tensor in tensor_list:
        mv = MultiVector(algebra, {0: tensor})  # Grade 0 scalar
        result.append(mv)
    return result



```

```python


to_numpy(x_flat[0])
```

```python


```

```python
x.cpu()
```

```python
# get first value of E 
#E[0,0,0,0,:]
bm(E[0][0][0][0])
```

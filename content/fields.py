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

def make_fields(F,
                speed = 50,
                vector_scale =2,
                e12_bounds = 6,
                e0_bounds  = 4,
                e3_bounds  = 0,
                e12_N      = 21,
                e0_N       = 100,
                e3_N       = 1,
                ):
 

    x = make_x(e12_bounds = e12_bounds,
            e0_bounds  = e0_bounds,
            e3_bounds  = e3_bounds,
            e12_N      = e12_N,
            e0_N       = e0_N,
            e3_N       = e3_N,
            )
    Fx = F(x) # make F-field directly
    frames     = render_frames(x, Fx, vector_scale=vector_scale)
    graph_func = make_graph_func(frames, speed)
 

    ## uncomment to animate
    return  graph_func 

 
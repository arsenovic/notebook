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
written 01/22/2026<br>
<a href="mailto:alex@810lab.com">alex@810lab.com</a><br>
</div>


# Operator Class

## Summary 
A Operator class is created.  An example application to the differential is explored . 

## Context 
Geometric calcus requires operators, which are functions that take and return functions.  Perhaps we can apply it to the 'differential' or 'a-derivative' 


For example :
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
 lets look at an implementation 


## Implementation

```python
class Operator:
    def __init__(self, right, left, *args, **kw):
        """
        A non-commutative operator that can be applied to functions

        Summary of operations
        - self * func: applies operator on the right
        - func * self: applies operator on the left
        - self ^ func: 1/2 * (right - left)
        - self | func: 1/2 * (right + left)

    
        
        Args:
            right: A function that takes a function f and returns the operator applied to f on the right
            left: A function that takes a function f and returns the operator applied to f on the left
            *args, **kw: Additional arguments to be passed to the right and left functions
        """
        self._right = right
        self._left = left

        self.right = lambda f: right(f, *args, **kw)
        self.left  = lambda f:  left(f, *args, **kw)


    def __call__(self, *args, **kw):
        """return a new Operator with the same left and right but updated args/kw"""
        return type(self)(self._right, self._left,*args, **kw)
    
    def __add__(self, op):
        """operator * func: applies operator on the right"""
        raise NotImplementedError("Addition of Operators is not implemented yet.")


    def __mul__(self, f):
        """operator * func: applies operator on the right"""
        return self.right(f)
    
    def __rmul__(self, f):
        """func * operator: applies operator on the left"""
        return self.left(f)
    
    def __xor__(self, f):
        """operator ^ func: antisymmetric combination = 1/2 * (self * f - f * self)"""
        return lambda *args, **kwargs: 0.5 * ((self*f)(*args, **kwargs) - (f*self)(*args, **kwargs))
    
    def __rxor__(self, f):
        """operator ^ func: antisymmetric combination = 1/2 * (f*self - self * f)"""
        return -self.__xor__(f)  # subtraction is antisymmetric

    
    def __or__(self, f):
        """operator | func: symmetric combination = 1/2 * (self * f + f * self)"""
        return lambda *args, **kwargs: 0.5 * ((self*f)(*args, **kwargs) + (f*self)(*args, **kwargs))
    
    def __ror__(self, f):
        """operator | func: symmetric combination = 1/2 * (f*self + self * f)"""
        return self.__or__(f)  # addition is symmetric

```

## Test 

```python
from kingdon import calculus as calc 
from kingdon import Algebra
sta = Algebra(1,3,start_index=0)
locals().update(sta.blades)

 

da = Operator(right = lambda f,a: calc.da(f, a, direction='forward'),
              left  = lambda f,a: calc.da(f, a, direction='backward'))


f = lambda x: 2*x**20   # function with high curvature to see differences in derivative approximations
x = sta.random_vector() # location at which to evaluate the derivative
a = sta.random_vector() # direction for the derivative

((da(a=a)|f)(x),  # d_a*f +f*d_a
calc.da(f=f, a=a, direction='central')(x),
calc.da(f=f, a=a, direction='forward')(x))
```

```python
((da(a=a)*f)*da(a=a))(x)
```

```python

```

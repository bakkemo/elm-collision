# Efficient Collision detection in Elm
This package implements the Gilbert-Johnson-Keerthi (GJK) collision detection algorithm for 2D convex objects.
It is quite efficient, usually converging within one or two iterations.

## Basic Usage

```haskell
import Collision exposing (..)

-- this is what polySupport looked like in 0.14 code
polySupport : List Pt -> Pt -> Pt
polySupport list d =
    let
        dotList = List.map (dot d) list
        decorated = (List.map2 (,)) dotList list
        (m, p) = List.maximum decorated -- maximum now returns a Maybe b
    in
        p


poly1 = [(-15,-10),(0,15),(12,-5)]
poly2 = [(-9,13),(6,13),(-2,22)]
collision 10 (poly1, polySupport) (poly2, polySupport) == True
````
**Note:** the first parameter is max recursion depth. It can easily be elided by defining an auxiliary helper
like `myCollision = collision 100`. Control over recursion depth can be useful when defining your own support
functions.

## API

```haskell
type alias Pt = (Float, Float)
type alias Mink a = (a, a -> Pt -> Pt)

collision : Int -> Mink a -> Mink b -> Bool
```
**Note:** a `Mink b` is a pair of: a boundary object of type `b`, and a suppport function of type
`f: b -> Pt -> Pt` which given a boundary object, and a direction vector (given by a Pt), produces
a point on the boundary furthest in the direction of the vector.

**example**
```haskell
polySupport [(-15,-10),(0,15),(12,5)] (1,0) == (12,5)
polySupport [(-15,-10),(0,15),(12,5)] (0,-1) == (-15,10)
```

You can define your own boundary objects and corresponding support functions, perhaps to handle
circles. Look in GJK.elm in bakkemo/umwelt for just such an example. It doesn't make sense for this
library to prescribe a boundary representation (for circles, or any OTHER object type).

**N.B.**

Determining if a point is inside an object is just a special case of this: (pt, (\a b -> a)) : Mink Pt is a
perfectly valid Minkowski object.

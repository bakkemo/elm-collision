# Efficient Collision detection in Elm
This package implements the Gilbert-Johnson-Keerthi (GJK) collision detection algorithm for 2D convex objects.
It is quite efficient, usually converging within one or two iterations.

This is a fork of
[bakkemo/elm-collision](https://github.com/bakkemo/elm-collision) by [Michael
Bakkemo](https://github.com/bakkemo). The original package does not use Maybe's
and feels kind of broken.

## Basic Usage

```elm
import Collision exposing (..)

-- Example polySupport function
dot : Pt -> Pt -> Float
dot (x1,y1) (x2,y2) = (x1*x2) + (y1*y2)

polySupport : List Pt -> Pt -> Maybe Pt
polySupport list d =
    let
        dotList = List.map (dot d) list
        decorated = (List.map2 (,)) dotList list
        max = List.maximum decorated
    in
        case max of
        Just (m, p) -> Just p
        _ -> Nothing
        
poly1 = [(-15,-10),(0,15),(12,-5)]
poly2 = [(-9,13),(6,13),(-2,22)]
collision 10 (poly1, polySupport) (poly2, polySupport) == Just True
````
**Note:** the first parameter to collision is max recursion depth. It can easily be elided by defining an auxiliary helper
like `myCollision = collision 100`. Control over recursion depth can be useful when defining your own support
functions.

## API

```elm
type alias Pt = (Float, Float)
type alias Mink a = (a, (a -> Pt -> Maybe Pt))

collision : Int -> Mink a -> Mink b -> Maybe Bool
```
**Note:** a `Mink b` is a pair of: a boundary object of type `b`, and a suppport function of type
`f: b -> Pt -> Maybe Pt` which given a boundary object, and a direction vector (given by a Pt), produces
a point on the boundary furthest in the direction of the vector.

**example**
```elm
polySupport [(-15,-10),(0,15),(12,5)] (1,0) == Just (12,5)
polySupport [(-15,-10),(0,15),(12,5)] (0,-1) == Just (-15,10)
```

You can define your own boundary objects and corresponding support functions, perhaps to handle
circles. Look in GJK.elm in bakkemo/umwelt for just such an example. It doesn't make sense for this
library to prescribe a boundary representation (for circles, or any OTHER object type).

**N.B.**

Determining if a point is inside an object is just a special case of this: (pt, (\a b -> a)) : Mink Pt is a
perfectly valid Minkowski object.

## Credits

 * [Michael Bakkemo](https://github.com/bakkemo)
 * [Firas Zaidan](https://github.com/zaidan)

## License

See `LICENSE` file.

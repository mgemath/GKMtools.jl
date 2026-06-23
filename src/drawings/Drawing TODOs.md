# TODO for drawings:

Add notion of a positive drawing: an admissible drawing where for each edge, the vector representing it in the drawing is a _positive_ multiple of its axial function value.

Positivity, convexity, and strong convexity are all properties of admissible drawings, but positivitity and convexity are independent of each other.

## TODO:

- ~~In line with the existing functions, add functions to check for existence positive drawings.~~
  Done: `positive_drawing_representative(G; strong=false)` in `positive_drawings.jl` returns a generic
  representative of the (unique) all-positive chamber together with its admissibility/convexity analysis,
  or `nothing` when no positive drawing exists.
- ~~Is it true that among the admissible drawing representatives there is at most one positive one?~~
  Yes: positivity pins down the single all-`+1` sign vector (up to the global `C`/`-C` identification),
  hence at most one chamber. Confirmed empirically in `test/drawings/test_positive_drawings.jl`.
![](m2m.png)

###### A general many-to-many relation manager with:

* Automatic serialization ($<<$ and $>>$) for all
* Automatic comparison $==$, $<$, $>$, $\leq$, $\geq$ and set-theory operations $\cap$, $\cup$ and $\setminus$ for all
  relations
* Classical symbolic operations for sparse addressing ($\mathtt{Transpose}$, $+$, $-$, $\times$)
* Generalized dependence lookup by node sets
* Permutation of generalized nodes and elements
* Compression of nodes and elements
* Canonical representation of symmetry groups
* Depth-first (DFS) search from nodes and elements
* DFS search for types
* Topological ordering of nodes and elements
* Topological ordering of types
* Lexicographical ordering of all dependence relations
* Automatic instantiation in a heterogeneous container by means of builders

  And much more

The idea is for the user to define a domain-specific "zoo", then create the fauna, and her own domain-specific functions
and avoid housekeeping.

It is written in C-style templated C++ by Pedro Areias (IST)

All thanks to AI and some colleagues


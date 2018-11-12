# Homology

This is the pure Python module for calculating homology groups of free chain complexes.It allows to calculate simplicial homology and cohomology groups. Also it can be used for calculating Khovanov homology of knots, Carter cohomology groups for 2-cocycle invariants of knots and so on.

## How to use

1. Import the module
```
import homology
```

2. Create chain complex object
```
cc = ChainComplex()
```

3. Add groups ot the complex
```
cc.add_group(i, dim)
```

4. Set boundary homeomorphisms
```
cc.add_boundary_map(i, matrix)
```

5. Get homology group
```
orders = cc.get_homology(i)
```

the result is a tuple (a1, a2, ...) and it means that H_i = Z_a1 + Z_a2 + ...

6. Get coordinates of the group summand generators
```
gen = cc.get_generator(i, summand)
```

the result is a vector in i-th chain group in default basis (this basis used for setting boundary maps).

7. For calculating cohomoly groups, convert the defind chain complex to the cochain complex
```
ccc = cc.generate_cochain_complex()
```

and calculate groups by using get_homology(i) method.

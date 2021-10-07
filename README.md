# Homology

This is the pure Python module for calculating homology groups of free chain complexes. It allows to calculate simplicial homology and cohomology groups. Also it can be used for calculating Khovanov homology of knots, Carter cohomology groups for 2-cocycle invariants of knots and so on.

## How to use

1. Import the module
```
import homology
```

2. Create the chain complex object
```
cc = ChainComplex()
```

3. Add groups to the complex
```
cc.add_group(i, dim)
```

4. Set boundary homomorphisms
```
cc.add_boundary_map(i, matrix)
```

5. Get homology group
```
orders = cc.get_homology(i)
```

the result is a tuple (a1, a2, ...) and it means that H_i = Z_a1 + Z_a2 + ... Order 0 means that this summand is infinite cyclic group.

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

## Compiled versions of matrix_helper module

By default homology calculations use pure Python version of functions from matrix_helper.py. But it's also possible to use compiled versions of it. MatrixHelper_Nim.pyd is a result of porting functions to nim language. MatrixHelper_Cpp.pyd is a port to c++. Port to c++ use MPIR library for calculations with big integers. There are only two functions, which use big integers: `reduce` and `get_matrix_kernel`. All other functions use int type of values.

There are no reasons to use compiled versions of the module. c++ version slightly faster, near x1.5 - x3 times. Nim version of `reduce` and `get_matrix_kernel` functions are slowly.

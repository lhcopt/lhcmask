# Tests performed on pymask

For pysixtrack we are using the branch: "feature/cc_from_mad_and_sixtrack"
Sixtracktools is on the master

We test the commit 89ef222ad8371704bb055185358ff133a42be6de

## Test 1 - b1 with bb legacy macros

We select:
```
mode = 'b1_with_bb_legacy_macros'
```
We exeute the python
```bash
python 004_pymask.py
```

We execute the reference:
```bash
 cd ../hl_lhc_collision
 python ../../unmask.py main.mask parameters_for_unmask.txt --run
```

We check that the output is identical
```bash
diff fc.2 ../hl_lhc_collision/fc.2
diff fc.3 ../hl_lhc_collision/fc.3
diff fc.8 ../hl_lhc_collision/fc.8
diff fc.16 ../hl_lhc_collision/fc.16
diff fc.34 ../hl_lhc_collision/fc.34
```

## Test 2 - b1 with new bb tools
We select:
```
mode = 'b1_with_bb'
```

We exeute the python
```bash
python 004_pymask.py
```

We execute the reference (can be reused from the previous test):
```bash
 cd ../hl_lhc_collision
 python ../../unmask.py main.mask parameters_for_unmask.txt --run
```

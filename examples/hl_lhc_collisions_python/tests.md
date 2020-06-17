# Tests performed on pymask

For pysixtrack we are using the branch: "feature/cc_from_mad_and_sixtrack"
Sixtracktools is on the master

## Test 1
```
mode = 'b1_with_bb_legacy_macros'
```

```bash
python 004_pymask.py
```

We execute the reference:
```bash
 cd ../hl_lhc_collision
 python ../../unmask.py main.mask parameters_for_unmask.txt --run
```

# Tests performed on pymask

For pysixtrack we are using the branch: "feature/cc_from_mad_and_sixtrack"
Sixtracktools is on the master

We test the commit 89ef222ad8371704bb055185358ff133a42be6de

## Test 1 - b1 with bb legacy macros

We execute the reference:
```bash
 cd ../hl_lhc_collision
 python ../../unmask.py main.mask parameters_for_unmask.txt --run
```

We setup the python version:
```bash
cd ../hl_lhc_collisions_python

```
We select:
```python
mode = 'b4_without_bb'
```
We run the python version:
```
python 004_pymask.py
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
```python
mode = 'b1_with_bb'
```

We exeute the python version:
```bash
python 004_pymask.py
```

We execute the reference (can be reused from the previous test):
```bash
 cd ../hl_lhc_collision
 python ../../unmask.py main.mask parameters_for_unmask.txt --run
```

Check using pysixtrack lines:
```bash
python t003_fc_to_fort.py
python t004_compare_pysixtrack_lines.py
```

**Check crab cavities:**

In ```t005_check_crabs.py``` set:
```python
# B1 ip5
beam_track = 'b1'
ip_choice = 5
plane = 'y'
phi_weak = Phi
phi_c_weak = Phi_c
```

Run check:
```bash
python t005_check_crabs.py
```
Minor deviations on crab orbit are observed due to errors and tuning.

Repeat for IP1.

## Test 3 - b4 without bb

We execute the reference:
```bash
 cd ../hl_lhc_collision_nobb_b4
 python ../../unmask.py main.mask parameters_for_unmask.txt --run
```

We setup the python version:
```bash
cd ../hl_lhc_collisions_python

```
We select:
```python
mode = 'b4_without_bb'
```
We run the python version:
```
python 004_pymask.py
```

We check that the output is identical
```bash
diff fc.2 ../hl_lhc_collision/fc.2
diff fc.3 ../hl_lhc_collision/fc.3
diff fc.8 ../hl_lhc_collision/fc.8
diff fc.16 ../hl_lhc_collision/fc.16
diff fc.34 ../hl_lhc_collision/fc.34
```

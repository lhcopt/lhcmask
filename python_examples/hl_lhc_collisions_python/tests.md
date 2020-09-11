# Tests performed on pymask

Sixtracktools is on the master
We test the commit e301c179d937843f15eb4bc01e647e15833fdc09

N.B. For the tests to be successful make sure that cpymad and madx correspond to the same (updated!) version.



## Test 1 - b1 with bb legacy macros

We execute the reference:
```bash
 cd ../../examples/hl_lhc_collision
 python ../../unmask.py main.mask parameters_for_unmask.txt --run
```

We setup the python version:
```bash
cd ../../python_examples/hl_lhc_collisions_python

```
We select:
```python
mode = 'b1_with_bb_legacy_macros'
```
We run the python version:
```
python 000_pymask.py
```

We check that the output is identical
```bash
diff fc.2 ../../examples/hl_lhc_collision/fc.2
diff fc.3 ../../examples/hl_lhc_collision/fc.3
diff fc.8 ../../examples/hl_lhc_collision/fc.8
diff fc.16 ../../examples/hl_lhc_collision/fc.16
diff fc.34 ../../examples/hl_lhc_collision/fc.34
```


## Test 2 - b1 with new bb tools
We select:
```python
mode = 'b1_with_bb'
```

We exeute the python version:
```bash
python 000_pymask.py
```

We execute the reference (can be reused from the previous test):
```bash
 cd ../../examples/hl_lhc_collision
 python ../../unmask.py main.mask parameters_for_unmask.txt --run
```

Configure ```t004_compare_pysixtrack_lines.py``` for the test:
```python
# Test b1
path_test = './'
type_test = 'sixtrack'
path_ref = '../hl_lhc_collision'
type_ref = 'sixtrack
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
python t003_fc_to_fort.py
python t005_check_crabs.py
```
Minor deviations on crab orbit are observed due to errors and tuning.

Repeat for IP1.

## Test 3 - b4 without bb

**WARNING!!!! In this mode the leveling mask module is not working, values are forced in python and the module is skipped**

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
diff fc.2 ../../examples/hl_lhc_collision_nobb_b4/fc.2
diff fc.3 ../../examples/hl_lhc_collision_nobb_b4/fc.3
diff fc.8 ../../examples/hl_lhc_collision_nobb_b4/fc.8
diff fc.16 ../../examples/hl_lhc_collision_nobb_b4/fc.16
diff fc.34 ../../examples/hl_lhc_collision_nobb_b4/fc.34
```

## Test 4 - b4 from b2 without bb

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
mode = 'b4_from_b2_without_bb'
```
We run the python version:
```
python 000_pymask.py
```

We check that the output is identical
```bash
diff fc.2 ../../examples/hl_lhc_collision_nobb_b4/fc.2
diff fc.3 ../../examples/hl_lhc_collision_nobb_b4/fc.3
diff fc.8 ../../examples/hl_lhc_collision_nobb_b4/fc.8
diff fc.16 ../../examples/hl_lhc_collision_nobb_b4/fc.16
diff fc.34 ../../examples/hl_lhc_collision_nobb_b4/fc.34
```

## Test 5 - b4 from b2 with bb
We setup the python version:
```bash
cd ../hl_lhc_collisions_python

```
We select:
```python
mode = 'b4_from_b2_with_bb'
```
We run the python version:
```
python 004_pymask.py
```

**As we don't have a reference we can only perform the check on the crab cavities:**

In ```t005_check_crabs.py``` set:
```python
# B4 ip1
beam_track = 'b4'
ip_choice = 1
plane = 'x'
phi_weak = -Phi
phi_c_weak = -Phi_c
```

Run check:
```bash
python t003_fc_to_fort.py
python t005_check_crabs.py
```
Minor deviations on crab orbit are observed due to errors and tuning.

Repeat for IP5 using:
```python
# B4 ip5
beam_track = 'b4'
ip_choice = 5
plane = 'y'
phi_weak = Phi
phi_c_weak = Phi_c
```

## To check a mad test againt its reference

```bash
diff fc.2 reference/fc.2
diff fc.3 reference/fc.3
diff fc.8 reference/fc.8
diff fc.16 reference/fc.16
diff fc.34 reference/fc.34
```


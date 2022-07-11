# Tests performed on pymask

N.B. For the tests to be successful make sure that cpymad and madx correspond to the same (updated!) version.


## Test 1 - b1 with bb legacy macros

We execute the reference:
```bash
 cd ../../examples/hl_lhc_collision
 python ../../unmask.py main.mask parameters_for_unmask.txt --run-cpymad
```

We setup the python version:
```bash
cd ../../python_examples/hl_lhc_collisions_python

```
In config.yaml we select:
```yaml
mode : b1_with_bb_legacy_macros
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
In config.yaml we select:
```yaml
mode : b1_with_bb
```

We exeute the python version:
```bash
python 000_pymask.py
```

We execute the reference (can be reused from the previous test):
```bash
 cd ../../examples/hl_lhc_collision
 python ../../unmask.py main.mask parameters_for_unmask.txt --run-cpymad
```
In the folder ```checks_and_doc```:
Configure ```t004_check_output_consistency.py``` for the test:
```python
# Test b1
tests = [
    {
        'test_name': 'B1 - pymask sixtrack input vs madx mask',
        'path_test': '../',
        'type_test': 'sixtrack',
        'path_ref': '../../../examples/hl_lhc_collision',
        'type_ref': 'sixtrack',
        'rtol': 3e-5,
        'atol': 1e-12,
        'strict': False,
    },
    {
        'test_name': 'B1 - pymask xsuite vs pymask sixtrack input',
        'path_test': '../xsuite/line_bb_dipole_not_cancelled.json',
        'type_test': 'xsuite',
        'path_ref': '../',
        'type_ref': 'sixtrack',
        'rtol': 4e-7,
        'atol': 1e-100,
        'strict': True,
    }
]
```

Check using xsuite lines:
```bash
python t003_fc_to_fort.py
python t004_check_output_consistency.py 
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
 cd ../../examples/hl_lhc_collision_nobb_b4
 python ../../unmask.py main.mask parameters_for_unmask.txt --run-cpymad
```

We setup the python version:
```bash
cd ../../python_examples/hl_lhc_collisions_python

```
In config.yaml we select:
```yaml
mode: b4_without_bb,
force_leveling:
  on_sep8: -0.03425547139366354
  on_sep2: 0.14471680504084292
```
(the separations are forced in order to be consitent with the mad-x test case for which the same values are impored, as the legacy leveling macro does not work in this case).

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

## Test 4 - b4 from b2 without bb

We execute the reference:
```bash
 cd ../../examples/hl_lhc_collision_nobb_b4
 python ../../unmask.py main.mask parameters_for_unmask.txt --run-cpymad
```

We setup the python version:
```bash
cd ../../python_examples/hl_lhc_collisions_python

```
In config.yaml we select:
```yaml
mode: b4_from_b2_without_bb
force_leveling: 
  on_sep8: -0.03425547139366354
  on_sep2: 0.14471680504084292
```
(the separations are forced in order to be consitent with the mad-x test case for which the same values are impored, as the legacy leveling macro does not work in this case).


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
cd ../../python_examples/hl_lhc_collisions_python

```
In config.yaml we select:
```python
mode: b4_from_b2_with_bb,
```
We run the python version:
```
python 000_pymask.py
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

## Test 7 and 8 - Check the matching and the tracking

Select on config.yaml

```yaml
beam: b1_without_bb
```
and run 

```bash 
python t003_fc_to_fort.py
python t007_check_orbit_and_lin_normal_form.py
python t008_check_against_sixtrack.py
```

then repeat for 

```yaml
beam: b1_with_bb
```

and

```yaml
beam: b4_from_b2_without_bb
```

and

```yaml
beam: b4_from_b2_with_bb
```



## To check a mad test againt its reference

```bash
diff fc.2 reference/fc.2
diff fc.3 reference/fc.3
diff fc.8 reference/fc.8
diff fc.16 reference/fc.16
diff fc.34 reference/fc.34
```


# Tests performed on pymask

N.B. For the tests to be successful make sure that cpymad and madx correspond to the same (updated!) version.


## Test 7 and 8 - Check the matching and the tracking

Select on config.yaml

```yaml
beam: b1_with_bb
```
and run 

```bash 
python t003_fc_to_fort.py
python t007_check_orbit_and_lin_normal_form.py
python t008_check_against_sixtrack.py
```

then repeat for 

```yaml
beam: b4_from_b2_with_bb
```

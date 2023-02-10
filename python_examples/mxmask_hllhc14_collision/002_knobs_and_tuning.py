import yaml
import json

import xtrack as xt

# Read config file
with open('config_knobs_and_tuning.yaml','r') as fid:
    configuration = yaml.safe_load(fid)

# Load collider
with open('collider_01_bb_off.json', 'r') as fid:
    collider = xt.Multiline.from_dict(json.load(fid))

# Load orbit correction configuration
with open('corr_co.json', 'r') as fid:
    co_corr_config = json.load(fid)

# Set all knobs (crossing angles, dispersion correction, rf, crab cavities,
# experimental magnets, etc.)
for kk, vv in configuration['knob_settings'].items():
    collider.vars[kk] = vv

# Build trackers
collider.build_trackers()

# Twiss before correction
twb1_before = collider['lhcb1'].twiss()
twb2_before = collider['lhcb2'].twiss(reverse=True)

# Tunings
for line_name in ['lhcb1', 'lhcb2']:
    knob_names = configuration['knob_names'][line_name]

    # Correct closed orbit
    print(f'Correcting closed orbit for {line_name}')
    collider[line_name].correct_closed_orbit(
                            reference=collider[line_name+'_co_ref'],
                            correction_config=co_corr_config[line_name])

    # Match coupling
    print(f'Matching coupling for {line_name}')
    collider[line_name].match(
        vary=[
            xt.Vary(name=knob_names['c_minus_knob_1'],
                    limits=[-0.5e-2, 0.5e-2], step=1e-5),
            xt.Vary(name=knob_names['c_minus_knob_2'],
                    limits=[-0.5e-2, 0.5e-2], step=1e-5)],
        targets=[xt.Target('c_minus', 0, tol=1e-4)])

    # Match tune and chromaticity
    print(f'Matching tune and chromaticity for {line_name}')
    collider[line_name].match(
        vary=[
            xt.Vary(knob_names['q_knob_1'], step=1e-8),
            xt.Vary(knob_names['q_knob_2'], step=1e-8),
            xt.Vary(knob_names['dq_knob_1'], step=1e-8),
            xt.Vary(knob_names['dq_knob_2'], step=1e-8),
        ],
        targets = [
            xt.Target('qx', configuration['qx'], tol=1e-4),
            xt.Target('qy', configuration['qy'], tol=1e-4),
            xt.Target('dqx', configuration['dqx'], tol=0.05),
            xt.Target('dqy', configuration['dqy'], tol=0.05)])




# Checks

# RF

# Crabs
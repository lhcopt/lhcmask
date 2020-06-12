import numpy as np
import bb_tools as bbt

install_lenses_in_sequence = bbt.install_lenses_in_sequence

def generate_bb_dataframes(mad,
    ip_names=['ip1', 'ip2', 'ip5', 'ip8'],
    numberOfLRPerIRSide=[25, 20, 25, 20],
    harmonic_number=35640,
    bunch_spacing_buckets=10,
    numberOfHOSlices=11,
    bunch_population_ppb=None,
    sigmaz_m=None,
    remove_dummy_lenses=True):

    for pp in ['circ', 'npart', 'gamma']:
        assert mad.sequence.lhcb1.beam[pp] == mad.sequence.lhcb2.beam[pp]

    circumference = mad.sequence.lhcb1.beam.circ
    madx_reference_bunch_charge = mad.sequence.lhcb1.beam.npart
    relativistic_gamma = mad.sequence.lhcb1.beam.gamma
    relativistic_beta = np.sqrt(1 - 1.0 / relativistic_gamma ** 2)
    if bunch_population_ppb is not None:
        bunch_charge_ppb = bunch_population_ppb
    else:
        bunch_charge_ppb = madx_reference_bunch_charge

    if sigmaz_m  is not None:
        sigt = sigmaz_m
    else:
        sigt = mad.sequence.lhcb1.beam.sigt

    bb_df_b1 = bbt.generate_set_of_bb_encounters_1beam(
        circumference, harmonic_number,
        bunch_spacing_buckets,
        numberOfHOSlices, bunch_charge_ppb, sigt,
        relativistic_beta, ip_names, numberOfLRPerIRSide,
        beam_name = 'b1',
        other_beam_name = 'b2')

    bb_df_b2 = bbt.generate_set_of_bb_encounters_1beam(
        circumference, harmonic_number,
        bunch_spacing_buckets,
        numberOfHOSlices, bunch_charge_ppb, sigt,
        relativistic_beta, ip_names, numberOfLRPerIRSide,
        beam_name = 'b2',
        other_beam_name = 'b1')

    # Generate mad info
    bbt.generate_mad_bb_info(bb_df_b1, mode='dummy')
    bbt.generate_mad_bb_info(bb_df_b2, mode='dummy')

    # Install dummy bb lenses in mad sequences
    bbt.install_lenses_in_sequence(mad, bb_df=bb_df_b1, sequence_name='lhcb1')
    bbt.install_lenses_in_sequence(mad, bb_df=bb_df_b2, sequence_name='lhcb2')

    # Use mad survey and twiss to get geometry and locations of all encounters
    bbt.get_geometry_and_optics_b1_b2(mad, bb_df_b1, bb_df_b2)

    # Get the position of the IPs in the surveys of the two beams
    ip_position_df = bbt.get_survey_ip_position_b1_b2(mad, ip_names)

    # Get geometry and optics at the partner encounter
    bbt.get_partner_corrected_position_and_optics(
            bb_df_b1, bb_df_b2, ip_position_df)

    # Compute separation, crossing plane rotation and crossing angle
    for bb_df in [bb_df_b1, bb_df_b2]:
        bbt.compute_separations(bb_df)
        bbt.compute_dpx_dpy(bb_df)
        bbt.compute_local_crossing_angle_and_plane(bb_df)

    # Get bb dataframe and mad model (with dummy bb) for beam 3 and 4
    bb_df_b3 = bbt.get_counter_rotating(bb_df_b1)
    bb_df_b4 = bbt.get_counter_rotating(bb_df_b2)
    bbt.generate_mad_bb_info(bb_df_b3, mode='dummy')
    bbt.generate_mad_bb_info(bb_df_b4, mode='dummy')

    # Generate mad info
    bbt.generate_mad_bb_info(bb_df_b1, mode='from_dataframe',
            madx_reference_bunch_charge=madx_reference_bunch_charge)
    bbt.generate_mad_bb_info(bb_df_b2, mode='from_dataframe',
            madx_reference_bunch_charge=madx_reference_bunch_charge)
    bbt.generate_mad_bb_info(bb_df_b3, mode='from_dataframe',
            madx_reference_bunch_charge=madx_reference_bunch_charge)
    bbt.generate_mad_bb_info(bb_df_b4, mode='from_dataframe',
            madx_reference_bunch_charge=madx_reference_bunch_charge)

    bb_dfs = {
        'b1': bb_df_b1,
        'b2': bb_df_b2,
        'b3': bb_df_b3,
        'b4': bb_df_b4}

    if remove_dummy_lenses:
        for beam in ['b1', 'b2']:
            bbdf = bb_dfs[beam]
            mad.input(f'seqedit, sequence={"lhc"+beam};')
            mad.input('flatten;')
            for nn in bbdf.elementName.values:
                print(f'remove, element={nn}')
                mad.input(f'remove, element={nn}')
            mad.input('flatten;')
            mad.input(f'endedit;')

    return bb_dfs

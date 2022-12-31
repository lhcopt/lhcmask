import xfields as xf
import pymask as pm

def install_dummy_bb_lenses(bb_df, line):

    ip_names = bb_df['ip_name'].unique().tolist()

    s_ips = {}
    for iipp in ip_names:
        s_ips[iipp] = line.get_s_position(iipp)

    for nn in bb_df.index:
        print(f'Insert: {nn}     ', end='\r', flush=True)
        ll = bb_df.loc[nn, 'label']
        iipp = bb_df.loc[nn, 'ip_name']

        if ll == 'bb_ho':
            new_bb = xf.BeamBeamBiGaussian3D(phi=0, alpha=0, other_beam_q0=0.,
                slices_other_beam_num_particles=[0],
                slices_other_beam_zeta_center=[0],
                slices_other_beam_Sigma_11=[1],
                slices_other_beam_Sigma_12=[0],
                slices_other_beam_Sigma_22=[0],
                slices_other_beam_Sigma_33=[1],
                slices_other_beam_Sigma_34=[0],
                slices_other_beam_Sigma_44=[0],
                )
        elif ll == 'bb_lr':
            new_bb = xf.BeamBeamBiGaussian2D(
                other_beam_beta0=1.,
                other_beam_q0=0,
                other_beam_num_particles=0.,
                other_beam_Sigma_11=1,
                other_beam_Sigma_33=1,
            )
        else:
            raise ValueError('Unknown label')

        line.insert_element(element=new_bb,
                                    at_s=(s_ips[bb_df.loc[nn, 'ip_name']]
                                        + bb_df.loc[nn, 'atPosition']),
                                    name=nn)

def install_beambeam_elements_in_lines(line_b1, line_b4, ip_names,
            circumference, harmonic_number, bunch_spacing_buckets,
            num_long_range_elems_per_side, num_slices_head_on,
            bunch_num_particles, sigmaz_m):

    # TODO: use keyword arguments
    # TODO: what happens if bunch length is different for the two beams
    bb_df_b1 = pm.generate_set_of_bb_encounters_1beam(
        circumference, harmonic_number,
        bunch_spacing_buckets,
        num_slices_head_on,
        bunch_num_particles, line_b1.particle_ref.q0,
        sigmaz_m, line_b1.particle_ref.beta0[0], ip_names, num_long_range_elems_per_side,
        beam_name = 'b1',
        other_beam_name = 'b2')


    bb_df_b2 = pm.generate_set_of_bb_encounters_1beam(
        circumference, harmonic_number,
        bunch_spacing_buckets,
        num_slices_head_on,
        bunch_num_particles, line_b4.particle_ref.q0,
        sigmaz_m,
        line_b4.particle_ref.beta0[0], ip_names, num_long_range_elems_per_side,
        beam_name = 'b2',
        other_beam_name = 'b1')
    bb_df_b2['atPosition'] = -bb_df_b2['atPosition'] # I am installing in b4 not in b2

    install_dummy_bb_lenses(bb_df=bb_df_b1, line=line_b1)
    install_dummy_bb_lenses(bb_df=bb_df_b2, line=line_b4)

    return bb_df_b1, bb_df_b2

def configure_beam_beam_elements(bb_df_b1, bb_df_b2, tracker_b1, tracker_b4,
                                 nemitt_x, nemitt_y, crab_strong_beam, ip_names):
    twiss_b1 = tracker_b1.twiss()
    twiss_b4 = tracker_b4.twiss()
    twiss_b2 = twiss_b4.reverse()

    survey_b1 = tracker_b1.survey()
    survey_b2 = tracker_b4.survey(reverse=True)

    sigmas_b1 = twiss_b1.get_betatron_sigmas(nemitt_x=nemitt_x, nemitt_y=nemitt_y)
    sigmas_b2 = twiss_b2.get_betatron_sigmas(nemitt_x=nemitt_x, nemitt_y=nemitt_y)

    # Use survey and twiss to get geometry and locations of all encounters
    pm.get_geometry_and_optics_b1_b2(
        mad=None,
        bb_df_b1=bb_df_b1,
        bb_df_b2=bb_df_b2,
        xsuite_line_b1=tracker_b1.line,
        xsuite_line_b2=tracker_b4.line,
        xsuite_twiss_b1=twiss_b1,
        xsuite_twiss_b2=twiss_b2,
        xsuite_survey_b1=survey_b1,
        xsuite_survey_b2=survey_b2,
        xsuite_sigmas_b1=sigmas_b1,
        xsuite_sigmas_b2=sigmas_b2,
    )

    # Get the position of the IPs in the surveys of the two beams
    ip_position_df = pm.get_survey_ip_position_b1_b2(mad=None, ip_names=ip_names,
        xsuite_survey_b1=survey_b1, xsuite_survey_b2=survey_b2)

    # Get geometry and optics at the partner encounter
    pm.get_partner_corrected_position_and_optics(
            bb_df_b1, bb_df_b2, ip_position_df)

    # Compute separation, crossing plane rotation, crossing angle and xma
    import pdb; pdb.set_trace()
    for bb_df in [bb_df_b1, bb_df_b2]:
        pm.compute_separations(bb_df)
        pm.compute_dpx_dpy(bb_df)
        pm.compute_local_crossing_angle_and_plane(bb_df)
        pm.compute_xma_yma(bb_df)

    # Get bb dataframe and mad model (with dummy bb) for beam 3 and 4
    bb_df_b3 = pm.get_counter_rotating(bb_df_b1)
    bb_df_b4 = pm.get_counter_rotating(bb_df_b2)

    bb_dfs = {
        'b1': bb_df_b1,
        'b2': bb_df_b2,
        'b3': bb_df_b3,
        'b4': bb_df_b4}

    if crab_strong_beam:
        pm.crabbing_strong_beam_xsuite(bb_dfs,
            tracker_b1, tracker_b4)
    else:
        print('Crabbing of strong beam skipped!')

    pm.setup_beam_beam_in_line(tracker_b1.line, bb_df_b1, bb_coupling=False)
    pm.setup_beam_beam_in_line(tracker_b4.line, bb_df_b4, bb_coupling=False)

    xf.configure_orbit_dependent_parameters_for_bb(tracker=tracker_b1,
                        particle_on_co=twiss_b1.particle_on_co)
    xf.configure_orbit_dependent_parameters_for_bb(tracker=tracker_b4,
                        particle_on_co=twiss_b4.particle_on_co)

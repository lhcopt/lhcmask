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
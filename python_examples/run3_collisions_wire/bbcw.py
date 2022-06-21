import numpy as np
import pandas as pd
import inspect

import pymask as pm
import pymask.pymasktools as pmTools


#=====================================================
# Constants for WIRES
WIRE_MARKERS = {}
WIRE_MARKERS['b1'] = ['bbcwe.4l5.u.b1',
                      'bbcwi.4l5.u.b1',
                      'bbcwe.4l5.d.b1',
                      'bbcwi.4l5.d.b1',
                      'bbcwb.a4l1.u.b1',
                      'bbcwt.a4l1.u.b1',
                      'bbcwb.a4l1.d.b1',
                      'bbcwt.a4l1.d.b1']
WIRE_MARKERS['b2'] = ['bbcwe.4r5.d.b2',
                      'bbcwe.4r5.u.b2',
                      'bbcwi.4r5.u.b2',
                      'bbcwi.4r5.d.b2',
                      'bbcwb.a4r1.d.b2',
                      'bbcwt.a4r1.d.b2',
                      'bbcwb.a4r1.u.b2',
                      'bbcwt.a4r1.u.b2']

# INSTALL, ELEMENT=string, CLASS=string,AT=real, FROM={string|SELECTED};

WIRE_INSTALL_FALLBACK = {}
WIRE_INSTALL_FALLBACK['b1'] =[('bbcwe.4l5.u.b1' ,'-148.54+(-308-ip5ofs.b1)*ds','ip5'),
                              ('bbcwi.4l5.u.b1' ,'-148.54+(-308-ip5ofs.b1)*ds','ip5'),
                              ('bbcwe.4l5.d.b1' ,'-147.35+(-308-ip5ofs.b1)*ds','ip5'),
                              ('bbcwi.4l5.d.b1' ,'-147.35+(-308-ip5ofs.b1)*ds','ip5'), 
                              ('bbcwb.a4l1.u.b1','-146.54+(0-ip1ofs.b1)*ds'   ,'ip1'),
                              ('bbcwt.a4l1.u.b1','-146.54+(0-ip1ofs.b1)*ds'   ,'ip1'),
                              ('bbcwb.a4l1.d.b1','-145.35+(0-ip1ofs.b1)*ds'   ,'ip1'),
                              ('bbcwt.a4l1.d.b1','-145.35+(0-ip1ofs.b1)*ds'   ,'ip1')]

WIRE_INSTALL_FALLBACK['b2'] = [('bbcwb.a4r1.d.b2','145.35+(0-ip1ofs.b2)*ds'  ,'ip1'),
                               ('bbcwt.a4r1.d.b2','145.35+(0-ip1ofs.b2)*ds'  ,'ip1'),
                               ('bbcwb.a4r1.u.b2','146.54+(0-ip1ofs.b2)*ds'  ,'ip1'),
                               ('bbcwt.a4r1.u.b2','146.54+(0-ip1ofs.b2)*ds'  ,'ip1'),                              
                               ('bbcwe.4r5.u.b2' ,'148.54+(308-ip5ofs.b2)*ds','ip5'),
                               ('bbcwi.4r5.u.b2' ,'148.54+(308-ip5ofs.b2)*ds','ip5'),
                               ('bbcwe.4r5.d.b2' ,'147.35+(308-ip5ofs.b2)*ds','ip5'),
                               ('bbcwi.4r5.d.b2' ,'147.35+(308-ip5ofs.b2)*ds','ip5')]




QFF_RECIPE = {('b1','ip1') : """kq4.r1b1 := kq4.r1b1_0 + enable_QFF*(-1.516687604126627e-12) *({I}/{d_w}^2);
                                kq4.l1b1 := kq4.l1b1_0 + enable_QFF*(-5.4967990643973015e-12)*({I}/{d_w}^2);""",
              
              ('b1','ip5') : """kq4.r5b1 := kq4.r5b1_0 + enable_QFF*(1.3654027047949919e-12) *({I}/{d_w}^2);
                                kq4.l5b1 := kq4.l5b1_0 + enable_QFF*(5.466219213479711e-12)  *({I}/{d_w}^2);""",
              
              ('b2','ip1') : """kq4.r1b2 := kq4.r1b2_0 + enable_QFF*(-5.494756708945384e-12) *({I}/{d_w}^2);
                                kq4.l1b2 := kq4.l1b2_0 + enable_QFF*(-1.5168382251914348e-12)*({I}/{d_w}^2);""",
              
              ('b2','ip5') : """kq4.r5b2 := kq4.r5b2_0 + enable_QFF*(5.464166813485956e-12)  *({I}/{d_w}^2);
                                kq4.l5b2 := kq4.l5b2_0 + enable_QFF*(1.3656162643083467e-12) *({I}/{d_w}^2);""", 
             }

ALIGN_RECIPE = {'b1':  {'bbcwe': """{xma} :=      {x_b}+{rw} ; {yma} := {y_b}     ;""",
                        'bbcwi': """{xma} :=      {x_b}-{rw} ; {yma} := {y_b}     ;""",
                        'bbcwt': """{xma} :=      {x_b}      ; {yma} := {y_b}+{rw};""",
                        'bbcwb': """{xma} :=      {x_b}      ; {yma} := {y_b}-{rw};"""},
                'b2':  {'bbcwe': """{xma} := -1*({x_b})+{rw} ; {yma} := {y_b}     ;""",
                        'bbcwi': """{xma} := -1*({x_b})-{rw} ; {yma} := {y_b}     ;""",
                        'bbcwt': """{xma} := -1*({x_b})      ; {yma} := {y_b}+{rw};""",
                        'bbcwb': """{xma} := -1*({x_b})      ; {yma} := {y_b}-{rw};"""}}

#=====================================================




# Utilities
#=====================================================
def install_wires(mad,configuration,seq_name):
    
    # installing monitors if not in seq_name:
    if configuration['wires_at_fallback']:
        to_install = pd.DataFrame(WIRE_INSTALL_FALLBACK[seq_name[-2:]],columns=['element','at','from'])
        to_install.insert(0,'mode','install')
        to_install.insert(2,'class','monitor')
        
        pmTools.seqedit(mad,seq_name=seq_name,editing = to_install)
    
    
    WIRE_LENGTH = configuration['wires_L']
    
    # Extracting information for each wire
    wires = [tuple(wire.split('.')) for wire in WIRE_MARKERS[seq_name[-2:]]]

    
    # Defining knobs
    knobs=[]
    for w_type,loc,s,beam in wires:
        knobs.append(f"""
        {w_type}_current.{loc}.{beam} = 0;
        {w_type}_xma.{loc}.{s}.{beam} = 1;
        {w_type}_yma.{loc}.{s}.{beam} = 1;""")
    madInput = '\n'.join(knobs)
    
    # Definition of each wire
    template = """{name} : wire,  
                    current := {current}/2,
                    L        = 0,
                    L_phy    = {L},
                    L_int    = 2*{L}, 
                    xma     := {xma}, 
                    yma     := {yma};
                """
    wires_def = [ template.format(name    = f'{w_type}_wire.{loc}.{s}.{beam}',
                                  current = f'{w_type}_current.{loc}.{beam}',
                                  xma     = f'{w_type}_xma.{loc}.{s}.{beam}',
                                  yma     = f'{w_type}_yma.{loc}.{s}.{beam}',
                                  L       = WIRE_LENGTH) for w_type,loc,s,beam in wires]

    # DEFINITION FOR INSTALLATION
    madInput += ';\n'.join(wires_def) + ';\n'

    # SEQEDIT FOR INSTALLATION
    madInput += pmTools.seqedit(mad,
                                madInput = False,
                                seq_name = seq_name,
                                editing  = {'mode'      : 'replace',
                                            'element'   : WIRE_MARKERS[seq_name[-2:]],
                                            'by'        : [_def.split(':')[0] for _def in wires_def]})
    mad.input(madInput)
    
    # Creating single knobs for the wires in each IP:
    make_knobs(mad)
#============================================================    

def remove_duplicates(myList):
    return list(set(myList))

#============================================================
def make_knobs(mad):
    
    # Extracting wire information for both beams
    wires = [tuple(wire.split('.')) for wire in WIRE_MARKERS['b1'] + WIRE_MARKERS['b2']]
    
    
    # Defining single current knob for the wires in each IP and single rw knob
    entries_current = []
    entries_rw      = []
    new_knobs       = []
    for w_type,loc,s,beam in wires:  
            #Saving knobs for later loop
            new_knobs.append((f'bbcw_I_ip{loc[-1]}.{beam}',f'bbcw_rw_ip{loc[-1]}.{beam}'))

            # MADCALL += Current
            entries_current.append(f"""bbcw_I_ip{loc[-1]}.{beam} = 0; {w_type}_current.{loc}.{beam} := bbcw_I_ip{loc[-1]}.{beam};""")
            # MADCALL += rw
            entries_rw.append(f'bbcw_rw_ip{loc[-1]}.{beam} = 1;' + ALIGN_RECIPE[beam][w_type].format(  xma = f'{w_type}_xma.{loc}.{s}.{beam}',
                                                                                                        yma = f'{w_type}_yma.{loc}.{s}.{beam}',
                                                                                                        rw  = f'bbcw_rw_ip{loc[-1]}.{beam}',
                                                                                                        x_b = 0,y_b=0))

    madCall   = '\n'.join(remove_duplicates(entries_current) + ['',''])
    madCall  += '\n'.join(entries_rw)
    
    mad.input(madCall)
#============================================================
    
#============================================================
def make_QFF_links(mad,configuration):
    mad.input(f"enable_QFF = {int(configuration['enable_QFF'])};")
    
    # Extracting wire information for both beams
    wires = [tuple(wire.split('.')) for wire in WIRE_MARKERS['b1'] + WIRE_MARKERS['b2']]

    all_knobs = []
    for w_type,loc,s,beam in wires:  
            all_knobs.append((f'bbcw_I_ip{loc[-1]}.{beam}',f'bbcw_rw_ip{loc[-1]}.{beam}'))
    all_knobs = remove_duplicates(all_knobs)



    # Linking kq4 values to wires for QFF!
    entries = []
    for _I,_rw in all_knobs:
        info = _I.split('_')[-1]
        ip,beam = info.split('.')

        recipe = inspect.cleandoc(QFF_RECIPE[(beam,ip)].format(I = _I, d_w = _rw))
        initialize = [f'{kq4}_0 = {kq4};' for kq4 in [line.split(' :=')[0] for line in recipe.split('\n')]]

        entries+= initialize + recipe.split('\n') + ['']

    madCall = '\n'.join(entries)
    
    mad.input(madCall)
#============================================================

#============================================================
def align_wires(mad,seq_name,_twiss=None):
    # TWISS to get beam position @ wire location
    if _twiss is None:
        mad.use(seq_name)
        mad.twiss()
        _twiss = mad.table.twiss.dframe()
    

    # Setting x and y positions of the wires with regard to the beam   
    entries = []
    for wire,wire_twiss in _twiss.groupby('keyword').get_group('wire').iterrows():
        
        w_type,loc,s,beam = wire.split('.')
        entries.append(ALIGN_RECIPE[beam][w_type[:-5]].format(xma = f'{w_type[:-5]}_xma.{loc}.{s}.{beam}',
                                                        yma = f'{w_type[:-5]}_yma.{loc}.{s}.{beam}',
                                                        rw  = f'bbcw_rw_ip{loc[-1]}.{beam}',
                                                        x_b = wire_twiss.x,
                                                        y_b = wire_twiss.y))
        
    madCall  = '\n'.join(entries)
    mad.input(madCall)
#============================================================


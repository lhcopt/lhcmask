import pymask as pm
import os

def prepare_pytrain_input(mad, configuration):
    # for the moment we assume from https://gitlab.cern.ch/agorzaws/train/-/blob/master/opt_train_2015.manb#L78 
    # not to slice the HO
    
    # To be clarified:
    # 1. the reflect of B2
    # 2. the sectorpure effect
    # 3. Is this needed?
    #  beamx = beam%lhcb2->ex;
    #  beamy = beam%lhcb2->ey;
    #  sigz  = beam%lhcb1->sigt;
    #  sige  = beam%lhcb1->sige;
    # 4. clarify the octupoles
    
    mad.input('''
    
    !beamx = beam%lhcb2->ex;
    !beamy = beam%lhcb2->ey;
    !sigz  = beam%lhcb1->sigt;
    !sige  = beam%lhcb1->sige;
    
    INSTALL_SINGLE_BB_MARK(label,BIM,nn,where,origin) : macro = {
    !! Install a single beam-beam marker
      install,element=bbmk_labelBIM_nn,class=bbmarker,at=where,from=origin;
    };
    
    INSTALL_BB_MARK(BIM) : macro = {
    !! Install all ho and parasitic beam-beam marker within +/-170 m from IP1/2/5/8 for a given beam
      option,warn,info;
      bbmarker: marker;
      Linteraux:=170.;
      nparasitic:=Linteraux/b_h_dist;
      seqedit,sequence=lhcBIM;
      
      ! Head-on
      ! IR1
      where=1.e-9;exec INSTALL_SINGLE_BB_MARK(ho1,BIM,0,where,IP1);
      n=1;while ( n <= (nho_IR1/2)) {exec HOLOC(1,$n);
      exec INSTALL_SINGLE_BB_MARK(ho.L1,BIM,$n,wherem,IP1.L1);
      exec INSTALL_SINGLE_BB_MARK(ho.R1,BIM,$n,wherep,IP1);
      n=n+1;};
      
      ! IR2
      where=1.e-9;exec INSTALL_SINGLE_BB_MARK(ho2,BIM,0,where,IP2);
      n=1;while ( n <= (nho_IR2/2)) {exec HOLOC(2,$n);
      exec INSTALL_SINGLE_BB_MARK(ho.L2,BIM,$n,wherem,IP2);
      exec INSTALL_SINGLE_BB_MARK(ho.R2,BIM,$n,wherep,IP2);
      n=n+1;};
      
      ! IR5
      where=1.e-9;exec INSTALL_SINGLE_BB_MARK(ho5,BIM,0,where,IP5);
      n=1;while ( n <= (nho_IR5/2)) {exec HOLOC(5,$n);
      exec INSTALL_SINGLE_BB_MARK(ho.L5,BIM,$n,wherem,IP5);
      exec INSTALL_SINGLE_BB_MARK(ho.R5,BIM,$n,wherep,IP5);
      n=n+1;};
      
      ! IR8
      where=1.e-9;exec INSTALL_SINGLE_BB_MARK(ho8,BIM,0,where,IP8);
      n=1;while ( n <= (nho_IR8/2)) {exec HOLOC(8,$n);
      exec INSTALL_SINGLE_BB_MARK(ho.L8,BIM,$n,wherem,IP8);
      exec INSTALL_SINGLE_BB_MARK(ho.R8,BIM,$n,wherep,IP8);
      n=n+1;};
      
      n=1;while ( n <= nparasitic) {where=-n*b_h_dist;exec INSTALL_SINGLE_BB_MARK(par.L1,BIM,$n,where,IP1.L1);n=n+1;};
      n=1;while ( n <= nparasitic) {where= n*b_h_dist;exec INSTALL_SINGLE_BB_MARK(par.R1,BIM,$n,where,IP1)   ;n=n+1;};
      n=1;while ( n <= nparasitic) {where=-n*b_h_dist;exec INSTALL_SINGLE_BB_MARK(par.L2,BIM,$n,where,IP2)   ;n=n+1;};
      n=1;while ( n <= nparasitic) {where= n*b_h_dist;exec INSTALL_SINGLE_BB_MARK(par.R2,BIM,$n,where,IP2)   ;n=n+1;};
      n=1;while ( n <= nparasitic) {where=-n*b_h_dist;exec INSTALL_SINGLE_BB_MARK(par.L5,BIM,$n,where,IP5)   ;n=n+1;};
      n=1;while ( n <= nparasitic) {where= n*b_h_dist;exec INSTALL_SINGLE_BB_MARK(par.R5,BIM,$n,where,IP5)   ;n=n+1;};
      n=1;while ( n <= nparasitic) {where=-n*b_h_dist;exec INSTALL_SINGLE_BB_MARK(par.L8,BIM,$n,where,IP8)   ;n=n+1;};
      n=1;while ( n <= nparasitic) {where= n*b_h_dist;exec INSTALL_SINGLE_BB_MARK(par.R8,BIM,$n,where,IP8)   ;n=n+1;};
      endedit;
      option,-warn,-info;
  };
    
    b_t_dist = 25;
    !n_inside_D1 = 5;
    b_h_dist := LHCLENGTH/HRF400 * 10./2. * b_t_dist / 25.;
    
    nho_IR1 = 0;
    nho_IR2 = 0;
    nho_IR5 = 0;
    nho_IR8 = 0;

    !call, file="beambeam_macros/macro_bb.madx";                  ! macros for beam-beam

    exec, INSTALL_BB_MARK(b1);
    exec, INSTALL_BB_MARK(b2);

    use, sequence=lhcb2;
    seqedit,sequence=lhcb2;
    flatten;
    reflect;
    endedit;
    use, sequence=lhcb2;
    
    use, sequence=lhcb1;
    seqedit,sequence=lhcb1;
    flatten;
    endedit;
    ''')
    
    
    mad.use('lhcb1')   
    pm.match_tune_and_chromaticity(mad,
        q1=configuration['qx0'],
        q2=configuration['qy0'],
        dq1=configuration['chromaticity_x'],
        dq2=configuration['chromaticity_y'],
        tune_knob1_name=configuration['knob_names']['qknob_1']['lhcb1'],
        tune_knob2_name=configuration['knob_names']['qknob_2']['lhcb1'],
        chromaticity_knob1_name=configuration['knob_names']['chromknob_1']['lhcb1'],
        chromaticity_knob2_name=configuration['knob_names']['chromknob_2']['lhcb1'],
        sequence_name='lhcb1', skip_use=True)
     
        
    mad.use('lhcb2')   
    pm.match_tune_and_chromaticity(mad,
        q1=configuration['qx0'],
        q2=configuration['qy0'],
        dq1=configuration['chromaticity_x'],
        dq2=configuration['chromaticity_y'],
        tune_knob1_name=configuration['knob_names']['qknob_1']['lhcb2'],
        tune_knob2_name=configuration['knob_names']['qknob_2']['lhcb2'],
        chromaticity_knob1_name=configuration['knob_names']['chromknob_1']['lhcb2'],
        chromaticity_knob2_name=configuration['knob_names']['chromknob_2']['lhcb2'],
        sequence_name='lhcb2', skip_use=True)
              
    mad.input(''' 
    
    select,flag=sectormap,clear;
    select,flag=sectormap,range=#e;
    select,flag=sectormap,class=marker,pattern=^bbmk_par.*;
    select,flag=sectormap,class=marker,pattern=^bbmk_ho.*_0;
    
    select,flag=twiss,clear;
    select,flag=twiss,class=marker,pattern=^ip.*, column=name,s,x,betx,alfx,dx,y,bety,alfy,dy;
    select,flag=twiss,class=marker,pattern=^bbmk_par.*, column=name,s,x,betx,alfx,dx,y,bety,alfy,dy;
    select,flag=twiss,class=marker,pattern=^bbmk_ho.*_0, column=name,s,x,betx,alfx,dx,y,bety,alfy,dy;

    select,flag=survey,clear;
    select,flag=survey,class=marker,pattern=^ip.* ,column=name,x,y,z;
    select,flag=survey,class=marker,pattern=^bbmk_par.*, column=name,x,y,z;
    select,flag=survey,class=marker,pattern=^bbmk_ho.*_0, column=name,x,y,z;

    set, format="22.18e";
    use, sequence=lhcb1;
    twiss,sequence=lhcb1, sectormap, sectorpure=true, sectorfile=train.mapf, file=train.optf;
    survey,sequence=lhcb1,file=train.surf;
    use, sequence=lhcb2;
    twiss,sequence=lhcb2, sectormap, sectorpure=true, sectorfile=train.mapb, file=train.optb;
    survey,sequence=lhcb2,file=train.surb;
    ''')
    
    os.system('rm -rf train-output')
    os.system('mkdir train-output')
    os.system("cat train.mapf | sed -E 's/BBMK_PAR.([LR])([0-9])B[12]_([0-9]+)/MKIP\\2P\\1\\3/ig' | sed -E 's/BBMK_HO([0-9])B[12]_0/MKIP\\1/ig' | tail -n +9 > train-output/train.manf")
    os.system("cat train.mapb | sed -E 's/BBMK_PAR.([LR])([0-9])B[12]_([0-9]+)/MKIP\\2P\\1\\3/ig' | sed -E 's/BBMK_HO([0-9])B[12]_0/MKIP\\1/ig' | tail -n +9 > train-output/train.manb")
    os.system("cat train.optf | sed -E 's/BBMK_PAR.([LR])([0-9])B[12]_([0-9]+)/MKIP\\2P\\1\\3/ig' | sed -E 's/BBMK_HO([0-9])B[12]_0/MKIP\\1/ig' > train-output/train.optf")
    os.system("cat train.optb | sed -E 's/BBMK_PAR.([LR])([0-9])B[12]_([0-9]+)/MKIP\\2P\\1\\3/ig' | sed -E 's/BBMK_HO([0-9])B[12]_0/MKIP\\1/ig' > train-output/train.optb")
    os.system("cat train.surf | sed -E 's/BBMK_PAR.([LR])([0-9])B[12]_([0-9]+)/MKIP\\2P\\1\\3/ig' | sed -E 's/BBMK_HO([0-9])B[12]_0/MKIP\\1/ig' > train-output/train.surf")
    os.system("cat train.surb | sed -E 's/BBMK_PAR.([LR])([0-9])B[12]_([0-9]+)/MKIP\\2P\\1\\3/ig' | sed -E 's/BBMK_HO([0-9])B[12]_0/MKIP\\1/ig' > train-output/train.surb")
    os.system("rm train*b")
    os.system("rm train*f")

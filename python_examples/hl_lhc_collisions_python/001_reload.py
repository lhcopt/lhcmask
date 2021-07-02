import os, shutil

import numpy as np

from cpymad.madx import Madx

def norm(x):
    return np.sqrt(np.sum(np.array(x) ** 2))

def trim_last_zeros(a):
    zz=np.where(a!=0)[0]
    if len(zz)>0:
        return a[:zz[-1]+1]
    else:
        return a



def fc_to_fort():
    for iff in [2,8,16,34]:
        os.system(f"rm fort.{iff}")
        try:
            shutil.copy(f"fc.{iff}", f"fort.{iff}")
        except Exception as ex:
            print(ex)
            print(f"fc.{iff} not found!")

    with open("fort.3", "w") as fout:
        with open("checks_and_doc/fort_parts/fort_beginning.3", "r") as fid_fort3b:
            fout.write(fid_fort3b.read())
        with open("fc.3", "r") as fid_fc3:
            fout.write(fid_fc3.read())
        with open("checks_and_doc/fort_parts/fort_end.3", "r") as fid_fort3e:
            fout.write(fid_fort3e.read())



mad=Madx()

mad.call("final_seq.madx")
mad.use(sequence="lhcb1")
mad.twiss()
mad.readtable(file="final_errors.tfs", table="errtab")
mad.seterr(table="errtab")
mad.set(format=".15g")
mad.twiss()
mad.sixtrack(cavall=True,radius=0.017,max_mult_ord=18,mult_auto_off=True)
fc_to_fort()

import sixtracktools
import xline

sixin=sixtracktools.sixinput.SixInput(".")
linesix=xline.Line.from_sixinput(sixin)
linemad=xline.Line.from_madx_sequence(mad.sequence.lhcb1,apply_madx_errors=True)

#linemad=xline.Line.from_madx_sequence(mad.sequence.lhcb1)

for ll in (linesix, linemad):
    ll.remove_inactive_multipoles(inplace=True)
    ll.remove_zero_length_drifts(inplace=True)
    ll.merge_consecutive_drifts(inplace=True)
    ll.merge_consecutive_multipoles(inplace=True)

for lb,li in [("fmad",linemad),("fsix",linesix)]:
    with open(lb,"w") as fh:
        for el,en in zip(li.elements,li.element_names):
           fh.write(str(en)+"\n")
           fh.write(str(el)+"\n")

len(linesix.elements)
len(linemad.elements)


for ii, (ee_test, ee_six, nn_test, nn_six) in enumerate(
    zip(linemad.elements, linesix.elements, linemad.element_names, linesix.element_names)
):
    assert type(ee_test) == type(ee_six)
    dtest = ee_test.to_dict(keepextra=True)
    dsix = ee_six.to_dict(keepextra=True)
    for kk in dtest.keys():
        if dtest[kk] == dsix[kk]:
            continue
        try:
            if hasattr(dtest[kk],"__iter__"):
                atest=trim_last_zeros(np.array(dtest[kk]))
                asix=trim_last_zeros(np.array(dsix[kk]))
            else:
                atest=(np.array(dtest[kk]))
                asix=(np.array(dsix[kk]))
            diff_rel = norm(atest - asix) / norm(dtest[kk])
        except ZeroDivisionError:
            diff_rel = 100.0
        if diff_rel < 1e-9:
            continue
        diff_abs = norm(atest - asix)
        if diff_abs > 0 and kk!="length":
            print(f"{nn_test}[{kk}] - test:{dtest[kk]} six:{dsix[kk]}")
        if diff_abs < 1e-12:
            continue




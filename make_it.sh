wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b  -p ./miniconda -f
source miniconda/bin/activate
python -m pip install ipython numpy scipy pandas psutil xsuite matplotlib
git clone -b release/v0.1.0 git@github.com:xsuite/tree_maker.git
python -m pip install -e tree_maker
git clone git@github.com:lhcopt/lhcmask.git -b release/v1.3.3
python -m pip install -e lhcmask
python -m pip install sixtracktools
python -m pip install NAFFlib
git clone git@github.com:lhcopt/lhcerrors.git
git clone git@github.com:lhcopt/lhctoolkit.git
git clone git@github.com:lhcopt/hllhc15.git
git clone git@github.com:lhcopt/hllhc14.git
# you need gfortran
# `sudo yum install gcc-gfortran`
# then `cd beambeam_macros`
# and `gfortran headonslice.f -o headonslice`
git clone git@github.com:lhcopt/beambeam_macros.git
git clone $(whoami)@lxplus.cern.ch:/afs/cern.ch/eng/lhc/optics/runIII

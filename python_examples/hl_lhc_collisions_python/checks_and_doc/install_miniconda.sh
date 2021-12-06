wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b  -p ./miniconda -f
rm Miniconda3-latest-Linux-x86_64.sh
source miniconda/bin/activate
pip install ipython jupyterlab numpy scipy pandas awkward matplotlib pyarrow pyyaml cpymad 

git clone https://gitlab.cern.ch/mihostet/pytrain.git
pip install -e pytrain

git clone https://github.com/xsuite/xobjects
pip install -e xobjects

git clone https://github.com/xsuite/xline
pip install -e xline

git clone https://github.com/xsuite/xpart
pip install -e xpart

git clone https://github.com/xsuite/xtrack
pip install -e xtrack

git clone https://github.com/xsuite/xfields
pip install -e xfields

git clone https://github.com/sixtrack/sixtracktools
pip install -e sixtracktools

git clone https://github.com/PyCOMPLETE/FillingPatterns.git
pip install -e FillingPatterns

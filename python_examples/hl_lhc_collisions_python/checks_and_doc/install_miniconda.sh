wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh         
bash Miniconda3-latest-Linux-x86_64.sh -b  -p ./miniconda -f                       
source miniconda/bin/activate                                                      
pip install ipython jupyterlab numpy scipy pandas awkward matplotlib               
pip install pyarrow pyyaml pytest                                                 
pip install cpymad                                                                 
pip install xsuite                                                                 
git clone git@github.com:lhcopt/lhcerrors.git                                      
git clone git@github.com:lhcopt/lhctoolkit.git                                     
git clone git@github.com:lhcopt/lhcmask.git   
git clone git@github.com:lhcopt/beambeam_macros.git
cd beambeam_macros
gfortran headonslice.f -o  headonslice
cd ..
cd lhcmask                                                                         
git checkout release/v1.3.1                                                        
pip install -e .                                                                   
cd ../                                                                             
git clone https://github.com/PyCOMPLETE/FillingPatterns.git                        
pip install ./FillingPatterns                                                      
python -m pip install sixtracktools            
python -m pip install NAFFlib

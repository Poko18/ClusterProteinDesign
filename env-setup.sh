#!/usr/bin/bash

##########################
#install Anaconda
##########################

#Installers are at https://repo.anaconda.com/archive/

anaconda_pkgname="Anaconda3-2024.02-1-Linux-x86_64.sh"
anaconda_base_path="${HOME:-/opt}/anaconda3"
#Download installer script
wget https://repo.anaconda.com/archive/$anaconda_pkgname

#Enable execute permission
chmod +x $anaconda_pkgname

#Execute anaconda installer script with -b (batch mode) for no manual intervention
./$anaconda_pkgname -b  -p $anaconda_base_path
source $anaconda_base_path/etc/profile.d/conda.sh

###########################
#create conda environments
###########################

#Create pyrosetta environment
conda create -n pyrosetta -c https://conda.rosettacommons.org pyrosetta

#Create colabdesign environment
conda create -n colabdesign python=3.9 pip
conda activate colabdesign
pip install https://storage.googleapis.com/jax-releases/cuda11/jaxlib-0.4.14+cuda11.cudnn86-cp39-cp39-manylinux2014_x86_64.whl
conda install pytorch==1.8.0 torchvision==0.9.0 torchaudio==0.8.0 cudatoolkit=11.1 -c pytorch -c conda-forge
pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1
#Downgrade to match cluster version of cuda
pip install nvidia-cuda-nvcc-cu11==11.7.99
pip install nvidia-cuda-runtime-cu11==11.7.99
conda install pyrosetta
pip install jax==0.4.14
pip install scipy==1.11.2
pip install dm-haiku==0.0.10
pip install chex==0.1.82
pip install optax==0.1.7
pip install flax==0.7.0
#Create RFdiffusion environment
#...
git clone https://github.com/RosettaCommons/RFdiffusion
cd RFdiffusion
mkdir models && cd models
wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt

conda env create -f env/SE3nv.yml
conda activate SE3nv
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../.. # change into the root directory of the repository
pip install -e . # install the rfdiffusion module from the root of the repository
pip install torch==1.9.1+cu111 torchvision==0.10.1+cu111 torchaudio==0.9.1 -f https://download.pytorch.org/whl/torch_stable.html
 pip install scipy==1.10.0

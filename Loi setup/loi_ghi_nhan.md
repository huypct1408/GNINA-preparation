1. /home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/environment_gnina.yml

  ```
name: gnina_env

channels:
  - pytorch
  - nvidia
  - conda-forge
  - defaults

dependencies:
  - python=3.11

  # Deep learning
  - pytorch
  - torchvision
  - torchaudio
  - pytorch-cuda=12.1

  # Cheminformatics
  - rdkit=2024.03.1
  - openbabel=3.1.1

  # Structure repair
  - pdbfixer
  - mdanalysis

  # Runtime
  - boost>=1.81
  - eigen>=3.4
  - glog
  - protobuf>=3.20

  # Scientific stack
  - numpy>=1.21
  - scipy
  - pandas
  - scikit-learn

  # Tools
  - matplotlib
  - jupyter
  - ipython
  - tqdm
  - git

  - pip
  - pip:
      - git+https://github.com/forlilab/molscrub.git
      - py3Dmol>=2.0.0
      - pyquaternion>=0.9.0
      - psutil>=5.8.0
```
I created a new venv, but 
'/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/bin/gnina' --help
/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/bin/gnina: error while loading shared libraries: libcudnn.so.9: cannot open shared object file: No such file or directory
conda list | grep cudnn
cudnn                     9.10.1.4             haad7af6_0    conda-forge
libcudnn                  9.10.1.4             h7d33bf5_0    conda-forge
libcudnn-dev              9.10.1.4             h0fdc2d1_0    conda-forge
pytorch                   2.5.1           py3.11_cuda12.1_cudnn9.1.0_0    pytorch

I download https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2 from GNINA
How can i fix the errors

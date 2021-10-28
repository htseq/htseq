#!/bin/bash
if [ $OS_NAME == 'ubuntu-latest' ]; then
  echo "Installing deps for linux"
  #sudo apt-get install -y swig3.0

elif [ $OS_NAME == 'macos-latest' ]; then
  echo "Find out OSX version"
  osx_version=$(sw_vers -productVersion)
  echo "OSX version: $osx_version"
  osx_ver1=$(echo $osx_version | cut -d. -f1)
  osx_ver2=$(echo $osx_version | cut -d. -f2)
  if [ $osx_ver1 -lt 10 ] || [ $osx_ver2 -lt 14 ]; then
    echo "OSX version not for deployment: $osx_version"
  else
    echo "OSX version for deployment: $osx_version"
  fi

  echo "Installing deps for OSX"
  # Prepare to exit upon failure
  set -e
  CONDA_URL="https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
  wget -nv "${CONDA_URL}"
  bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda
  echo "$PATH"
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate

  # Make conda environment and activate
  conda create -y -n ci python=$CONDA_PY
  conda activate ci

  # Use pip from conda
  conda install -y pip
  pip --version

  # Use SWIG from conda
  conda install -c conda-forge swig

else
  echo "OS not recognized: $OS_NAME"
  exit 1
fi

echo "Install Python dependencies"
pip install setuptools wheel pytest Cython numpy pysam matplotlib

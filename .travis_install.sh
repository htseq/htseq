#!/bin/bash
# try and make wheels
if [ $DOCKER_IMAGE ]; then
  docker run --rm -v `pwd`:/io $DOCKER_IMAGE /io/buildwheels.sh 
  if [ $? != 0 ]; then
      exit 1
  fi
  ls wheelhouse/
  if [ $? != 0 ]; then
      exit 1
  fi

# compile normally
else
  if [ $TRAVIS_OS_NAME == 'osx' ]; then
    export PATH="$HOME/miniconda/bin:$PATH"
    source $HOME/miniconda/bin/activate
  fi
  
  # Requirements should be installed by pip
  #if [ $TRAVIS_OS_NAME == 'linux' ]; then
  #  sed -i "s|pysam>=0.9.0|$PYSAM_VERSION|" requirements.txt
  #elif [ $TRAVIS_OS_NAME == 'osx' ]; then
  #  sed -i "" "s|pysam>=0.9.0|$PYSAM_VERSION|" requirements.txt
  #else
  #  echo "OS not recognized: $TRAVIS_OS_NAME"
  #  exit 1
  #fi
  #if [ $? != 0 ]; then
  #    exit 1
  #fi

  echo "Installing requirements"
  pip install -r requirements.txt
  echo "Requirements installed"
  
  if [ -n "${PYPI}" ]; then
    echo "Installing HTSeq from testpypi"
    pip install -i "${PYPI}" HTSeq
  else
    echo "Installing HTSeq from production pypi"
    pip install HTSeq
  fi
  if [ $? != 0 ]; then
      exit 1
  fi
  echo "HTSeq installed from pypi"
fi

#!/bin/bash

PYTHON=python3

###############################################################################

set -e

CLEAN=0
CONDA=0
PYTEST_ARGS=test
VENV_DIR=.venv
VERBOSE=0
SKIP_INSTALL=0

while getopts ":coavst:k:" OPTION; do
    echo "$OPTION"
    case $OPTION in
        c)
          CLEAN=1
          ;;
	a)
	  CONDA=1
	  ;;
        o)
          PYTEST_ARGS=test/test_htseq-count.py
          ;;
        t)
	  PYTEST_ARGS=$OPTARG
          ;;
	k)
	  PYTEST_ARGS="${PYTEST_ARGS} -k $OPTARG"
	  ;;
        s)
          SKIP_INSTALL=1
          ;;
        v)
          VERBOSE=1
          ;;
        \?)
          echo "Usage: $0 [-coavtk]"
          ;;
    esac
done
shift $((OPTIND -1))


if [ x$CLEAN = x1 ]; then
    rm -rf build/
fi

if [ x$CONDA = x1 ]; then
  source /opt/anaconda/bin/activate
  conda activate scanpy
  PYTHON=python
  PIP=pip
  PYTEST=pytest
else
  PYTHON=$VENV_DIR/bin/python
  PIP=$VENV_DIR/bin/pip
  PYTEST=$VENV_DIR/bin/pytest

  if [ ! -d $VENV_DIR ]; then
      $PYTHON -m venv $VENV_DIR
  fi
  $VENV_DIR/bin/pip install -U pip wheel numpy pybigwig
fi

if [ x$SKIP_INSTALL = x0 ]; then
  $PIP install --use-feature=in-tree-build .[htseq-qa,test]
elif [ x$VERBOSE = x1 ]; then
  echo "Skipping install"
fi

if [ x$CONDA = x1 ]; then
  if [ x$VERBOSE = x1 ]; then
    echo "${PYTEST} ${PYTEST_ARGS}"
  fi
  $PYTEST ${PYTEST_ARGS}
else
  if [ x$VERBOSE = x1 ]; then
    echo "PATH=${VENV_DIR}/bin:${PATH} ${PYTEST} ${PYTEST_ARGS}"
  fi
  PATH=${VENV_DIR}/bin:${PATH} $PYTEST ${PYTEST_ARGS}
fi

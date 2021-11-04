#!/bin/bash
if [ $OS_NAME == 'macos-latest' ]; then
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate
  conda activate ci
fi

echo "Unit tests"
pytest test
if [ $? != 0 ]; then
    exit 1
fi

echo "Doctests on documentation"
pytest --doctest-glob="*.rst" doc/*.rst doc/tutorials/tss.rst
if [ $? != 0 ]; then
    exit 1
fi

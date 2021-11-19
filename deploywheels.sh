#!/bin/bash
#
# Build manylinux wheels for HTSeq. Based on the example at
# <https://github.com/pypa/python-manylinux-demo>
#
# It is best to run this in a fresh clone of the repository!
#
# Run this within the repository root:
#   docker run --rm -v $(pwd):/io quay.io/pypa/manylinux2010_x86_64 /io/buildwheels.sh
#
# The wheels will be put into the wheelhouse/ subdirectory.
#
# For interactive tests:
#   docker run -it -v $(pwd):/io quay.io/pypa/manylinux2010_x86_64 /bin/bash

set -xeuo pipefail

# only deploy builds for a release_<sematic-version>_RC?? tag to testpypi
echo "Figure out if a release is appropriate for this tag: ${GITHUB_REF}"
if [ -z $GITHUB_REF ]; then
  echo 'No GITHUB_REF, exit'
  exit 0
fi
TAG=$(echo $GITHUB_REF | cut -f3 -d/)
TAG1=$(echo $TAG | cut -f1 -d_)
TAG2=$(echo $TAG | cut -f2 -d_)
TAG3=$(echo $TAG | cut -f3 -d_)
if [ -z $TAG2 ]; then
  echo 'No TAG2, exit'
  exit 0;
fi
if [ $TAG1 != 'release' ] || [ $TAG2 != $(cat /io/VERSION) ]; then
  echo 'No release tag or wrong version, exit'
  exit 0;
fi

# deploy onto pypitest unless you have no RC
if [ -z $TAG3 ]; then
  TWINE_PASSWORD=${TWINE_PASSWORD_PYPI}
  TWINE_REPOSITORY='https://upload.pypi.org/legacy/'
  echo 'Deploying to production pypi'
elif [ ${TAG3:0:2} == 'RC' ]; then
  TWINE_PASSWORD=${TWINE_PASSWORD_TESTPYPI}
  TWINE_REPOSITORY='https://test.pypi.org/legacy/'
  echo 'Deploying to testpypi'
else
  echo "Tag not recognized: $GITHUB_REF"
  exit 1
fi

# Deploy binary packages
echo "PYTHON_VERSION: ${PYTHON_VERSION}"
export PYTHON_FDN=cp$(echo ${PYTHON_VERSION} | sed 's/\.//')
PYBINS="/opt/python/${PYTHON_FDN}*/bin"

# New manylinux spec means we have to query the actual line
HTSEQ_VERSION=$(cat /io/VERSION)
HTSEQ_WHEEL_FILE=$(ls /io/wheelhouse/HTSeq-*.whl)
echo "HTSEQ_VERSION: ${HTSEQ_VERSION}"
echo "HTSEQ_WHEEL_FILE: ${HTSEQ_WHEEL_FILE}"

ERRS=0
for PYBIN in ${PYBINS}; do
  PYVER=$(basename $(dirname ${PYBIN}))
  echo "PYVER=$PYVER"
  echo "TWINE_REPOSITORY=$TWINE_REPOSITORY"
  echo "TWINE_USERNAME=$TWINE_USERNAME"
  echo "TWINE_PASSWORD=$TWINE_PASSWORD"
  ${PYBIN}/pip install twine
  ${PYBIN}/twine upload --repository-url "${TWINE_REPOSITORY}" -u "${TWINE_USERNAME}" -p "${TWINE_PASSWORD}" ${HTSEQ_WHEEL_FILE}
  if [ $? != 0 ]; then
    ERRS=1
  fi
done

echo "Deploy source code only from one version"
if [ x$SOURCE_VERSION == x$PYTHON_VERSION ]; then
  ${PYBIN}/twine upload --repository-url "${TWINE_REPOSITORY}" -u "${TWINE_USERNAME}" -p "${TWINE_PASSWORD}" /io/wheelhouse/HTSeq-${HTSEQ_VERSION}.tar.gz
  if [ $? != 0 ]; then
    ERRS=1
  fi
fi
exit $ERRS

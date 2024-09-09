#!/bin/bash
if [ $OS_NAME == 'macos-latest' ]; then
	# only deploy builds for a release_<sematic-version>_RC?? tag to testpypi
	if [ -z $GITHUB_REF ]; then
		echo 'No GITHUB_REF, exit'
		exit 0
	fi
	TAG=$(echo $GITHUB_REF | cut -f3 -d/)
	TAG1=$(echo $TAG | cut -f1 -d_)
	TAG2=$(echo $TAG | cut -f2 -d_)
	TAG3=$(echo $TAG | cut -f3 -d_)
	echo "GITHUB_REF: ${GITHUB_REF}, TAG: ${TAG}, TAG1: ${TAG1}, TAG2: ${TAG2}, TAG3: ${TAG3}"
	if [ -z $TAG2 ]; then
		echo 'No TAG2, exit'
		exit 0
	fi
	if [ $TAG1 != 'release' ] || [ $TAG2 != $(cat VERSION) ]; then
		echo 'No release tag or wrong version, exit'
		exit 0
	fi

	# do not deploy on linux outside of manylinux
	if [ -z $DOCKER_IMAGE ] && [ $OS_NAME != 'macos-latest' ]; then
		echo 'Not inside manylinux docker image and not OSX, exit'
		exit 0
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

	echo "Deploying for OSX"
	# Prepare to exit upon failure
	set -e

	# Only deploy on 14.5 to ensure backwards compatibility
	osx_version=$(sw_vers -productVersion)
	echo "OSX version: $osx_version"
	osx_ver1=$(echo $osx_version | cut -d. -f1)
	osx_ver2=$(echo $osx_version | cut -d. -f2)
	if [ $osx_ver1 -lt 14 ] || [ $osx_ver2 -lt 5 ]; then
		echo "OSX version not for deployment (<14.5): $osx_version"
		exit 1
	fi

	HTSEQ_VERSION=$(cat VERSION)
	echo "TWINE_REPOSITORY=$TWINE_REPOSITORY"
	echo "TWINE_USERNAME=$TWINE_USERNAME"
	echo "TWINE_PASSWORD=$TWINE_PASSWORD"
	export PATH="$HOME/miniconda/bin:$PATH"
	source $HOME/miniconda/bin/activate
	conda activate ci

	# make wheel
	mkdir wheelhouse
	pip wheel . -w wheelhouse/
	if [ $? != 0 ]; then
		exit 1
	fi

	echo "Contents of wheelhouse:"
	ls wheelhouse

	echo "Figure out architecture string for wheel..."
	PYVER=$(echo $CONDA_PY | sed 's/\.//')
	PYARCH=cp${PYVER}-cp${PYVER}
	TWINE_WHEEL=$(ls wheelhouse/HTSeq-${HTSEQ_VERSION}-${PYARCH}*.whl)
	echo "TWINE_WHEEL=$TWINE_WHEEL"

	echo "Install twine for upload..."
	pip --version
	pip install twine

	echo "Uploading..."
	twine upload --repository-url "${TWINE_REPOSITORY}" -u "${TWINE_USERNAME}" -p "${TWINE_PASSWORD}" "${TWINE_WHEEL}"
	if [ $? != 0 ]; then
		echo "Upload of wheel failed"
		exit 1
	fi
	echo "Upload of wheel complete"

else
	echo "No DOCKER_IMAGE and not OSX, we should not be here!"
	exit 1
fi

# Analyzers

## Installation

Setup CMSSW environment and get Analyzer repository:

	cmsrel CMSSW_8_0_26_patch1
    cd CMSSW_8_0_26_patch1/src/
    cmsenv
    git cms-init

For MET Correction:
	
	git cms-merge-topic cms-met:METRecipe_8020 -u
	git cms-merge-topic cms-met:METRecipe_80X_part2 -u

Electron MVA ID:

	git cms-merge-topic ikrav:egm_id_80X_v2

Analysis package:

	git clone https://github.com/zctao/Analyzers.git

MiniAODHelper:

	git clone https://github.com/cms-ttH/MiniAOD.git
	(currently use branch CMSSW_8_0_24_v1_sync)

LeptonID package shared with ND ttH-Multilepton group:

	git clone https://github.com/cms-ttH/ttH-LeptonID.git ttH/LeptonID

Compile:

	scram b -j 16

Add the area containing the Electron MVA weights:

	cd $CMSSW_BASE/external
	cd slc6_amd64_gcc530/
	git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
	cd data/RecoEgamma/ElectronIdentification/data
	git checkout egm_id_80X_v1
	cd $CMSSW_BASE/src

When running with CRAB one needs to add the following option to the crab config file: config.JobType.sendExternalFolder = True This is needed until the PR including this ID will be integrated in CMSSW/cms-data.

## Formatting

The following code style(s) is(are) used:

* For C/C++: Linux Kernel Style, see
https://www.kernel.org/doc/Documentation/CodingStyle

* For Python: PEP, see https://www.python.org/dev/peps/pep-0008/

* **RUN clang-format**: some automatization for C++ is available through
clang-format. Use it! Run:

	$ clang-format -style=file -i *.h

	$ clang-format -style=file -i *.cc

	$ clang-format -style=file -i *.cpp

to update all files in a current directory.
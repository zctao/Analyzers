# Analyzers

## Installation

Setup CMSSW environment and get Analyzer repository:

	 cmsrel CMSSW_8_0_25
     cd CMSSW_8_0_25/src/
     cmsenv
     git cms-init

     git clone https://github.com/zctao/Analyzers.git

For MET Correction:
	
	git cms-merge-topic cms-met:METRecipe_8020

Analysis package:

	git clone https://github.com/zctao/Analyzers.git

MiniAODHelper:

	git clone https://github.com/cms-ttH/MiniAOD.git
	(currently use branch CMSSW_8_0_8)

LeptonID package shared with ND ttH-Multilepton group:

	git clone https://github.com/cms-ttH/ttH-LeptonID.git ttH/LeptonID

Fix needed for a possible bug in current version of CMSSW:

	git cms-addpkg CommonTools/Utils
	sed -i 's|Math/include|Math/interface|' CommonTools/Utils/interface/normalizedPhi.h

Compile:

	scram b -j 16

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
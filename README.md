# Analyzers

## Installation

Setup CMSSW environment and get Analyzer repository:

	 cmsrel CMSSW_7_6_3_patch2
     cd CMSSW_7_6_3_patch2/src/
     cmsenv
     git cms-init

     git clone https://github.com/zctao/Analyzers.git
	 cd Analyzers
   	 git checkout cmssw_7_6_X
   	 cd ..

Get dependencies:

	 git clone https://github.com/cms-ttH/MiniAOD.git

Then either:

	 git checkout CMSSW_7_6_3
	 (compatibility with LeptonID package is not guaranteed)

or stay in the master branch, but need to delete plugins (not used for now anyway) in both MiniAODHelper and SkimDilep directories due to a bug

For ttH, H->tautau, get LeptonID package shared with ttH-Multilepton group:
	git clone https://github.com/cms-ttH/ttH-LeptonID.git

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
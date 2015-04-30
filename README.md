# Analyzers

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
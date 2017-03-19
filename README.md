# About


TinkerTests is currently a standalone test suite for [Tinker](https://dasher.wustl.edu/tinker/), designed with the hope that it will find its way into the [main Tinker code](https://github.com/jayponder/tinker) eventually.

## Usage

This wrapper is designed to be very lightweight and easy to use, relying primarily on deposited input files and reliable reference outputs.  Full usage details, including parallel running and running subsets of tests, can be obtained by running

	runtests.py --help
	
from any directory.  To run all tests, simply call the `runtests.py` script from any directory, after setting the `TINKERDIR` environmental variable to a folder containing the Tinker binaries.

## Input markup

The script relies on special markers in the input file to determine how the calculation is to be run.  See the current tests for examples of how to do this.

#### #Description:
This should contain a (80 or fewer character) description of what this test does, for printing purposes.

#### #Labels:
One or more labels used to categorize this test case, to allow easy running of a subset of tests using the -L 
flag.

#### #Tolerance:
Specifies the precision to which certain quantities are checked.  Currently supported values:-

* Energy    (default is 1E-6)
* Gradient  (default is 1E-6)
* Geometry  (default is 1E-6)

These are specified as, *e.g.* `#Precision Energy 1E-6` with as many lines as necessary to fully specify all tolerances.

#### #RunCommand:
The command to run that will execute the test case.

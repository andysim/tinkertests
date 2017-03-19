#!/usr/bin/env python

from glob import glob
import os
import re
import subprocess
from multiprocessing.dummy import Pool as ThreadPool
import argparse
import sys

default_tolerances = {
    'energy'   : 1e-6,
    'gradient' : 1e-6,
    'geometry' : 1e-6,
}

try:
    TINKERDIR = os.environ['TINKERDIR']
except KeyError:
    TINKERDIR = None

#
# Parse options
#
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--numprocs',
                    type=int,
                    help='The number of tests to run simultaneously.',
                    default=1)
parser.add_argument('-l', '--labels',
                    nargs='*',
                    help='Ron only tests containing the label(s) specified')
parser.add_argument('-m', '--makerefs',
                    action='store_true',
                    help='Replaces reference files for tests (filtered by labels, if present) instead of testing.')
parser.add_argument('-f', '--file',
                    nargs=1,
                    help='Process only the key file specified in this command.')
parser.add_argument('-t', '--tinkerdir',
                    help='The location of the tinker binaries (overrides the default found in $TINKERDIR).',
                    default=TINKERDIR)
parser.add_argument('-v', '--verbose',
                    action='store_true',
                    help='Provide detailed output.')
args = parser.parse_args()

if not args.tinkerdir:
    raise Exception("Tinker executable directory not specied.  Set the TINKERDIR variable, or use the --tinkerdir flag.")

if not os.path.isdir(args.tinkerdir):
    raise Exception("Tinker executable directory not valid: %s doesn't exist." % args.tinkerdir)



CSI = "\x1B["
def make_green(s):
    return CSI + "32;1m" + s + CSI + "0m"

def make_red(s):
    return CSI + "31;1m" + s + CSI + "0m"


def run_tinker(commands, outfile, args):
    """ A helper utility to execute commands (specified as a list), writing the output
        to outfile, and returning the output as a list."""

    if not isinstance(commands, list):
        raise ValueError("The commands argument to run_tinker should be a list of strings")
    mycmd = commands[:]
    mycmd[0] = args.tinkerdir + "/" + mycmd[0]
    if args.verbose:
        print("Attempting to run: %s" % " ".join(mycmd))
    with open(outfile, 'w') as fp:
        process = subprocess.Popen(mycmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout,stderr = process.communicate()
        # Python 3 returns a byte stream, so we convert here, just in case
        stdout = stdout.decode(sys.stdout.encoding)
        stderr = stderr.decode(sys.stdout.encoding)
        fp.write(stdout)
        fp.write(stderr + "\n")
    if process.returncode:
        raise RuntimeError("Command\n\n\t%s\n\nFailed.  See %s for details." % (" ".join(commands), outfile))
    return stdout.splitlines()


def parse_keywords(keyfile):
    keywords = { 'tolerance'   : default_tolerances,
                 'description' : "",
                 'runcommand'  : "",
                 'labels'      : []}
    tolerances = default_tolerances.keys()
    for line in keyfile:
        line = line.strip().split()
        if not line: continue
        line[0] = line[0].lower()
        if line[0].startswith('#description'):
            keywords['description'] = " ".join(line[1:])
        elif line[0].startswith('#tolerance'):
            try:
                if len(line) != 3:
                    raise Exception()
                key = line[1]
                val = float(line[2])
                if key not in tolerances:
                    raise Exception()
                keywords['tolerance'][key] = val
            except:
                raise SyntaxError("Bad tolerance argument!  Should be\n\n#tolerance %s value\n\n" % "/".join(tolerances))
        elif line[0].startswith('#labels'):
            keywords['labels'] = line[1:]
        elif line[0].startswith('#runcommand'):
            keywords['runcommand'] = line[1:]

    if not keywords['runcommand']:
        raise SyntaxError('#RunCommand not specified.')
    return keywords


def validate(refvals, outvals, keywords, args):
    if refvals.keys() != outvals.keys():
        raise Exception("Different keys detected in outputs")
    failed = False
    for quantity in refvals.keys():
        if args.verbose:
            print("\tChecking %s" % quantity)
        try:
            tolerance = keywords['tolerance'][quantity]
        except KeyError:
            raise Exception("Comparison for %s is not supported." % quantity)
        if refvals[quantity].keys() != outvals[quantity].keys():
            raise Exception("Different keys detected for %s" % quantity)
        for component in refvals[quantity].keys():
            out = outvals[quantity][component]
            ref = refvals[quantity][component]
            this_failed = False
            for o,r in zip(out, ref):
                if abs(o-r) > tolerance:
                    failed = True
                    this_failed = True
            if args.verbose:
                if this_failed:
                    print(make_red("\t\t%s failed" % component))
                else:
                    print("\t\t%s passed" % component)
    return failed


def parse_testgrad(out):
    #Total Potential Energy :                 -0.31418519 Kcal/mole
    totpotre = re.compile(r' Total Potential Energy :\s*(-?\d+\.\d+) Kcal/mole')
    # Anlyt       1       0.40103733      0.43512325      0.35325000      0.68916525
    agradre = re.compile(r' Anlyt\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(?:-?\d+\.\d+)')
    # Numer       1       0.40103733      0.43512325      0.35325000      0.68916525
    ngradre = re.compile(r' Numer\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(?:-?\d+\.\d+)')
    # This is a tough one to compile, but we're going to use it to parse out stuff like
    # Potential Energy Breakdown by Individual Components :
    #
    #  Energy       EB              EA              EBA             EUB
    #  Terms        EAA             EOPB            EOPD            EID
    #               EIT             ET              EPT             EBT
    #               EAT             ETT             EV              EC
    #               ECD             ED              EM              EP
    #               ER              ES              ELF             EG
    #               EX
    #
    #           39.64034584    206.70961139      0.00000000     -4.31167395
    #            0.00000000      0.00000000      0.00000000      0.00000000
    #            0.00000000      0.00000000      0.00000000      0.00000000
    #            0.00000000      0.00000000     -0.16604375      0.00000000
    #            0.00000000      0.00000000     -0.21064571     -0.31418519
    #            0.00000000      0.00000000      0.00000000      0.00000000
    #            0.00000000
    #
    eheaderre = re.compile(r' Potential Energy Breakdown by Individual Components :\n?')
    elabellinere = re.compile(r'  (Energy|Terms)?(\s+E[A-Z]+)+\n?')
    evaluelinere = re.compile(r'(\s+(-?\d+\.\d+))+\n?')
    elabelre = re.compile(r'E[A-Z]+')
    evaluere = re.compile(r'-?\d+\.\d+')
    def parse_energy_terms(out):
        parsed_some_values = False
        labels = []
        values = []
        for line in out:
            # Look for labels
            match = elabellinere.match(line)
            if match:
                labels.extend(elabelre.findall(line))
                continue
            # Look for numbers
            match = evaluelinere.match(line)
            if match:
                values.extend(map(float, evaluere.findall(line)))
                parsed_some_values = True
                continue
            # The first blank line after the values signals the end for us
            if not line.strip() and parsed_some_values:
                return zip(labels, values)

    results = {}
    energy = {}
    gradient = {}
    agrad = []
    ngrad = []
    elabs = []
    evals = []
    for n,line in enumerate(out):
        # Total potential energy
        matches = totpotre.match(line)
        if matches: energy['Total Potential'] = [float(matches.group(1))]
        # Energy components
        matches = eheaderre.match(line)
        if matches:
            for l,v in parse_energy_terms(out[n:]):
                energy[l] = [v]
        # Analytic gradient
        matches = ngradre.match(line)
        if matches: ngrad.extend(map(float, matches.groups(1)))
        # Numerical gradient
        matches = agradre.match(line)
        if matches: agrad.extend(map(float, matches.groups(1)))

    if agrad: gradient['Analytic'] = agrad
    if ngrad: gradient['Numerical'] = ngrad
    results['energy'] = energy
    results['gradient'] = gradient
    return results


def parse_analyze(out):
    floatre = re.compile(r'-?\d+\.\d+')
    results = {}
    energy = {}
    gradient = {}
    geometry = {}
    for n,line in enumerate(out):
        sl = line.split()
        # dE/dV (Virial-based) :                  0.068212 Kcal/mole/A**3
        if line.startswith(" dE/dV (Virial-based) :"): gradient['Analytic Virial'] = [ float(sl[3]) ]
        # dE/dV (Finite Diff) :                   0.068208 Kcal/mole/A**3
        if line.startswith(" dE/dV (Finite Diff) :"): gradient['Numerical Virial'] = [ float(sl[4]) ]
        # Pressure (Analytical, 0 K) :           -4677.200 Atmospheres
        if line.startswith(" Pressure (Analytical, 0 K) :"): gradient['Analytical Pressure'] = [ float(sl[5]) ]
        # Pressure (Numerical, 0 K) :            -4676.892 Atmospheres
        if line.startswith(" Pressure (Numerical, 0 K) :"): gradient['Numerical Pressure'] = [ float(sl[5]) ]
        # Bond Stretching                        6349.1257          16569
        if line.startswith(" Bond Stretching                  "): energy['Bond Stretching'] = [ float(sl[2]) ]
        # Angle Bending                          3749.0542          11584
        if line.startswith(" Angle Bending                    "): energy['Angle Bending'] = [ float(sl[2]) ]
        # Stretch-Bend                            -19.7897           4031
        if line.startswith(" Stretch-Bend                     "): energy['Stretch-Bend'] = [ float(sl[1]) ]
        # Urey-Bradley                            668.1058           7023
        if line.startswith(" Urey-Bradley                     "): energy['Urey-Bradley'] = [ float(sl[1]) ]
        # Out-of-Plane Bend                       110.8436           1566
        if line.startswith(" Out-of-Plane Bend                "): energy['Out-of-Plane Bend'] = [ float(sl[2]) ]
        # Torsional Angle                         433.6959           6701
        if line.startswith(" Torsional Angle                  "): energy['Torsional Angle'] = [ float(sl[2]) ]
        # Pi-Orbital Torsion                       58.7507            292
        if line.startswith(" Pi-Orbital Torsion               "): energy['Pi-Orbital Torsion'] = [ float(sl[2]) ]
        # Torsion-Torsion                         -41.9998            147
        if line.startswith(" Torsion-Torsion                  "): energy['Torsion-Torsion'] = [ float(sl[1]) ]
        # Van der Waals                         31548.0594        8304545
        if line.startswith(" Van der Waals                    "): energy['Van der Waals'] = [ float(sl[3]) ]
        # Atomic Multipoles                    -78760.2884        1667721
        if line.startswith(" Atomic Multipoles                "): energy['Atomic Multipoles'] = [ float(sl[2]) ]
        # Polarization                         -31796.4657        1667721
        if line.startswith(" Polarization                     "): energy['Polarization'] = [ float(sl[1]) ]
        # Internal Virial Tensor :               17258.896     115.886    -307.770
        #                                          115.886   16300.180     778.428
        #                                         -307.770     778.428   15756.319
        if line.startswith(" Internal Virial Tensor :"):
            gradient['Virial Tensor'] = map(float, floatre.findall("".join(out[n:n+3])))
        # Dipole Moment Magnitude :               1244.807 Debyes
        if line.startswith(" Dipole Moment Magnitude :"):  energy['Dipole Norm'] = [ float(sl[4]) ]
        # Dipole X,Y,Z-Components :               1123.325     247.506     475.843
        if line.startswith(" Dipole X,Y,Z-Components :"): energy['Dipole Components'] = map(float, sl[3:6])
        # Quadrupole Moment Tensor :             -1402.777    3261.667    3912.700
        #      (Buckinghams)                      3261.667     900.860   12633.678
        #                                         3912.700   12633.678     501.917
        if line.startswith(" Quadrupole Moment Tensor :"):
            energy['Quadrupole Components'] = map(float, floatre.findall("".join(out[n:n+3])))
        #  Radius of Gyration :                      30.915 Angstroms
        if line.startswith(" Radius of Gyration :"): geometry['Radius of Gyration'] = [ float(sl[4]) ]
        # Total Potential Energy :             -67700.9082 Kcal/mole
        if line.startswith(" Total Potential Energy :"): energy['Total Potential'] = [ float(sl[4]) ]
    results['energy'] = energy
    results['gradient'] = gradient
    results['geometry'] = geometry
    return results


def check_results(command, out, ref, keywords, args):
    """  The general strategy here is to parse out results into a common data structure, a dict of dict of lists:
                results[type][name] = [ val(s) ]
         where type is energy, gradient, coords, etc. (see keys of default_tolerances above); name is
         the specific quantity being tests, e.g. polarization energy, bonded energy.  This way we can
         write a very generic validation routine (energies, forces and gradients can be tested the same way)
         that can obey tolerances provided by the user.  Note that this means energies are provided as a list.
         By storing the names in the second index, we can print out meaningful messages with verbose on. """

    if command[0] == 'testgrad':
        if args.verbose: print("Checking testgrad outputs")
        refvals = parse_testgrad(ref)
        outvals = parse_testgrad(out)
    elif command[0] == 'analyze':
        if args.verbose: print("Checking analyze outputs")
        refvals = parse_analyze(ref)
        outvals = parse_analyze(out)
    else:
        raise Exception("No handler defined to check %s yet!" % command[0])
    return validate(refvals, outvals, keywords, args)


def run_testcase(testcase):
    """ Runs all steps needed for a single test case """
    reffile = testcase + '.ref'
    outfile = testcase + '.out'

    with open('%s.key'%testcase, 'r') as fp:
        keywords = parse_keywords(fp)

    if args.labels and testcase not in args.labels:
        if args.verbose:
            print("Specified label not found: moving on")

    if args.makerefs:
        # Just make the reference outputs; no comparison
        print("\tUpdating reference output for %s" % testcase)
        run_tinker(keywords['runcommand'], reffile, args)
    else:
        # Run tests and compare to reference outputs
        if args.verbose:
            print("Working on %s...\n" % testcase)

        output = run_tinker(keywords['runcommand'], outfile, args)
        if check_results(keywords['runcommand'], output, open(reffile).readlines(), keywords, args):
            line = ' {0:.<86}FAILED'.format(keywords['description'])
            print(make_red(line))
            failures.append(testcase)
        else:
            line = ' {0:.<86}PASSED'.format(keywords['description'])
            print(line)


#
# Run the tests
#
try:
    # Change directory to the tests
    scriptpath = os.path.dirname(os.path.realpath(__file__))
    testspath = scriptpath + '/tests'
    os.chdir(testspath)

    # Figure out the list of work
    if args.file:
        testcases = [ os.path.splitext(args.file[0])[0] ]
    else:
        testcases = [ os.path.splitext(testname)[0] for testname in glob('*.key') ]

    print("\n\tTesting binaries located in %s" % args.tinkerdir)
    print("\tRunning from %s, using %d cores\n" % (testspath, args.numprocs))

    failures = []
    if args.makerefs:
        print("\t###########################################################")
        print("\t#  WARNING! Reference files will be updated.  You should  #")
        print("\t#  only be doing this with a reliable version of Tinker.  #")
        print("\t###########################################################")


    # Set up a pool of workers to crank out the work in parallel
    pool=ThreadPool(args.numprocs)
    pool.map(run_testcase, sorted(testcases))
    pool.close()
    pool.join()
except KeyboardInterrupt:
    print("\nTesting interrupted...\n")
    if len(failed):
        print("\n The following tests failed:\n")
        print("\n".join(failures))
    sys.exit(1)

if len(failures):
    print("\nThe following %d of %d tests failed:\n" % (len(failures), len(testcases)))
    print("\n".join(failures))
else:
    if not args.makerefs:
        print("\n All %d tests succeeded!\n" % len(testcases))

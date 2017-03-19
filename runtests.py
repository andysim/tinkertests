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


def check_results(command, out, ref, keywords, args):
    if command[0] == 'testgrad':
        if args.verbose: print("Checking testgrad outputs")
        refvals = parse_testgrad(ref)
        outvals = parse_testgrad(out)
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

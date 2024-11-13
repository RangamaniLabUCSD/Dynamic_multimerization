#!/usr/bin/env python3
#
# compare.py
#
# copyright F. Nedelec, 14.12.2007, 14.03.2018; 4.8.2020, 18.10.2023

"""
compare.py
    Compare source files contained in two root directories
Usage:
    compare.py root1 root2 [opendiff]
"""

import sys, os, shutil, time, tempfile
import difflib, subprocess

mode="diff"
diff="diff --side-by-side -W200 -p --suppress-common-lines"



def print_spacer(arg):
    """print a line of width size, with 'arg' in the middle"""
    rows, cols = os.popen('stty size', 'r').read().split()
    sys.stdout.write(chr(27)+"[36;2m"); sys.stdout.flush()
    print(arg.center(int(cols), '-'))
    sys.stdout.write(chr(27)+"[0m"); sys.stdout.flush()


def compareFiles(fileL, fileR):
    if not os.path.isfile(fileL):
        return print("missing < %s" % fileL)
    if not os.path.isfile(fileR):
        return print("missing > %s" % fileR)
    #print("compare %s %s" % (fileL, fileR))
    comp = os.popen(diff+" "+fileL+" "+fileR)
    difLR = comp.read()
    comp.close()
    if difLR:
        if mode == 'diff':
            print_spacer(os.path.basename(fileL))
            for line in difLR:
                print(line, end='')
            
            sys.stdout.write(chr(27)+"[33;1m")
            print("This was `%s`" % os.path.basename(fileL), end=':')
            ans = input(' (return/left/right/open/q) ?>'+chr(27)+'[0m')
            
            if ans == "left" or ans == "l":
                shutil.copyfile(fileL, fileR)
            elif ans == "right" or ans == "r":
                shutil.copyfile(fileR, fileL)
            elif ans == "swap":
                fid, tmp = tempfile.mkstemp('.txt', 'temp', '', True)
                print(os.getcwd(), os.path.isfile(file))
                os.close(fid);
                os.rename(fileL, tmp)
                os.rename(fileR, fileL)
                os.rename(tmp, fileR)
            elif ans == "open":
                os.system("opendiff "+fileL+" "+fileR+"&")
            elif ans == "q":
                sys.exit()
        else:
            subprocess.call(["opendiff", fileL, fileR])
            #we wait a bit for the application to start
            time.sleep(0.5)


def interesting(file):
    return ( file.endswith('.py') or file.endswith('.m')
        or file.endswith('.h') or file.endswith('.cc')
        or file.endswith('.md') or file.endswith('.txt')
        or file.startswith('makefile') or file.startswith('.cym') )


def process_dir(roots, pathL, files):
    """compare files in the current directory"""
    if pathL.startswith(roots[0]):
        path = pathL[len(roots[0]):]
    else:
        path = pathL
    pathR = os.path.normpath(roots[1] + '/' + path)
    #print("PATH %s %s %s" % (path, pathL, pathR))
    #print("FILES %s" % files)
    if path.endswith('.svn'):
        return
    if path.endswith('.git'):
        return
    if 'DerivedData/' in path:
        return
    if 0 <= path.find('/.git/'):
        return
    if path.startswith('/bin'):
        return
    if path.startswith('/build/'):
        return
    print_spacer('%s %s'%(pathL, pathR))
    for file in files:
        if interesting(file):
            fileL = os.path.join(pathL, file)
            fileR = os.path.join(pathR, file)
            compareFiles(fileL, fileR)

#------------------------------------------------------------------------------

def main(args):
    """main"""
    global mode
    if len(args) < 2:
        print("Error: you must specify two root directories!")
        sys.exit()
    rootL = os.path.relpath(args[0])
    rootR = os.path.relpath(args[1])
    if not os.path.isdir(rootL):
        print("Error: `%s' is not a directory" % rootL)
        sys.exit()
    if not os.path.isdir(rootR):
        print("Error: `%s' is not a directory" % rootR)
        sys.exit()
    #parse command-line arguments:    
    for arg in args[2:]:
        if arg == 'opendiff':
            mode=arg
        else:
            print("unknown argument '%s'" % arg)
            sys.exit()
    # process all directories within 'root':
    #print("Comparing %s and %s" % (rootL, rootR))
    for path, dirs, files in os.walk(rootL, topdown=False):
        process_dir([rootL, rootR], path, files)



if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])


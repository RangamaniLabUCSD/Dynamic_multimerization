#!/usr/bin/env python3
#
# make_image.py:
# create the image for one frame in cytosim's result file
#
# copyright F. Nedelec, 2008-2010
#

"""
    Calls executable to generate an image in each specified directory

Syntax:

    make_image.py executable-with-arguments [options] [directories]
 
    The current directory is used if none is specified.
     
Options are specified as 'option=value', without space around the '=' sign.
Existing options and their values:

    format    png, ppm             file format (default = 'png')
    output    file-name            name of file (default = image.png)
    lazy      0 or 1 (default=1)   bail out if output file exists

Examples: 

    make_image.py 'play frame=10' run*
    make_image.py 'play frame=100 window_size=512,256' lazy=0 run*

F. Nedelec, March 19 2011 and before. Revised Sept-Nov 2012.
"""


try:
    import sys, os, subprocess
except ImportError:
    print("  Error: could not load necessary python modules\n")
    sys.exit()


# some parameters:
playexe = []
formats = [ 'png', 'ppm', 'tga' ]
format  = 'png'
output  = 'image.'+format
lazy    = 1

err = sys.stderr

#------------------------------------------------------------------------

def is_image(s):
    return s.startswith('image') and s.endswith(format)


def makeImages(format):
    """call executable to generate one image in the current directory"""
    val = subprocess.call(playexe + ['image', 'image_format='+format], stderr=None)
    if val != 0:
        err.write("`%s' failed with value %i\n" % (playexe[0], val))
        return [];
    files = []
    with os.scandir() as it:
        for e in it:
            if e.is_file() and is_image(e.name):
                files.append(e.name)
    if not files:
        raise IOError("'%s' did not produce any image!" % playexe[0])
    #newest = max(images, key=lambda x: os.stat(x).st_mtime)
    return files


def makeImagesUnzip(format, input_file='objects.cmo'):
    """unzip cytosim's input if necessary to make the images"""
    rmfile=False
    if not os.path.isfile(input_file):
        input_file_gz = input_file + '.gz'
        if os.path.isfile(input_file_gz):
            cmd = "gunzip %s -c > %s" % (input_file_gz, input_file)
            subprocess.call(cmd, shell=True)
            rmfile=True
    if not os.path.isfile(input_file):
        err.write("file '%s' not found!\n" % input_file)
        return []
    images = makeImages(format)
    if rmfile and os.path.isfile(input_file_gz):
        os.remove(input_file)
    return images

#------------------------------------------------------------------------

def process(dirpath, filenames):
    """make image in directory dirpath"""
    if ( output in filenames ) and lazy:
        print("Folder '%s' already has %s" % (dirpath, output))
        return
    try:
        image = makeImagesUnzip(format)
        if len(image) == 0:
            err.write("no image produced in `%s'\n"%dirpath)
        elif len(image) == 1:
            os.rename(image[0], output)
            print(image[0] + ' -> ' + output)
    except Exception as e:
        err.write("Error in `%s': %s\n" % (dirpath, repr(e)))


def process_dir(dirpath):
    """call process() with appropriate arguments"""
    files = []
    with os.scandir() as it:
        for e in it:
            if e.is_file():
                files.append(e.name)
    process(dirpath, files)


#------------------------------------------------------------------------

def main(args):
    """process command line arguments"""
    global playexe, format, output, lazy, codec, colors
    paths = []
    
    playexe = args[0].split()
    arg0 = os.path.expanduser(playexe[0])
    if os.access(arg0, os.X_OK) and os.path.isfile(arg0):
        playexe[0] = os.path.abspath(arg0)
    else:
        err.write("Error: you must specify an executable on the command line\n")
        sys.exit()
    
    for arg in args[1:]:
        [key, equal, value] = arg.partition('=')
        
        if key=='' or equal!='=' or value=='':
            if os.path.isdir(arg):
                paths.append(arg)
            else:
                err.write("ignored '%s' on command line\n" % arg)
        else:
            if key=='output':
                output = value
            elif key=='format':
                format = value
                output = 'image.'+format
                if not format in formats:
                    err.write("Invalid image format '%s'\n" % format)
            elif key=='colors':
                colors = int(value)
            elif key=='lazy':
                lazy = int(value)
            else:
                err.write("ignored '%s' on command line\n" % arg)
    
    if not paths:
        paths.append('.')
    
    cwd = os.getcwd()
    for p in paths:
        print("make_image.py "+"- "*32+p)
        os.chdir(p)
        process_dir(p)
        os.chdir(cwd)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])



#!/usr/bin/env python3
#
# make_movie.py creates MPEG4 movies for cytosim
#
# Copyright F. J. Nedelec, 2007 - 2023
#
# To make MP4 Quicktime movies, you need to install ffmpeg:
# http://www.ffmpeg.org
# via Macports: sudo port install ffmpeg
# via Brew: brew install ffmpeg
#
# Note that given a set of images, one can assemble them directly with ffmpeg
# Movie with a subjective quality (lower is better), use -crf X, with X in [0, 51]:
# ffmpeg -r 12 -i image%04d.png -pix_fmt yuv420p -c:v libx264 -crf 23 movie.mp4
#
# Movie with fixed bit rate:
# ffmpeg -r 12 -i image%04d.png -pix_fmt yuv420p -c:v libx264 -b:v 2000k movie.mp4
#

"""
    Create movies for cytosim in given directories

Syntax:

    make_movie.py executable-with-arguments [options] [directories]

Procedure:
    
    In each directory, this will:
      1- call the executable to generate images
      2- invoke ffmpeg or other tools to assemble a movie
      3- delete images created in step 1
    
    If images of the right format are already present in the directory,
    you may specify an empty executable and step 1 and 3 will be skipped.

    The current directory is used if none is specified.
    
Options are specified as 'option=value', without space around the '=' sign.
Existing options and their values:

    codec     mpeg4, h264, h265   movie codec (default = 'h264')
    rate      integer             images per second (default = 12)
    quality   integer (default=3) subjective quality for MPEG4: 1=great, 4=good
                                  data rate for Quicktime (eg. 64)
    
    cleanup   0 or 1 (default=1)  remove temporary files (default = 1)
    lazy      0 or 1 (default=1)  bail out if output file already exists

Examples:

    make_movie.py 'play window_size=512,256' run*
    make_movie.py '~/bin/play3 zoom=2' run*
    make_movie.py imgs codec=265

History:
    Created by F. Nedelec, 14.12.2007
    Improved by Beat Rupp, March 2010
    Revised on March 19 2011 and Sept-Nov 2012 by F. Nedelec.
    F. Nedelec, 10.2013: images are created in sub-directories
    F. Nedelec, 12.2013, 9.2014, 04.2016, 15.03.2020, 10.01.2022
    
    https://ffmpeg.org/ffmpeg.html
"""

try:
    import sys, os, shutil, subprocess
except ImportError:
    print("  Error: could not load necessary python modules\n")
    sys.exit()

def executable(arg):
    return os.path.isfile(arg) and os.access(arg, os.X_OK)

# some parameters:
tool    = []
codec   = 'h264'
rate    = 12
lazy    = 1
cleanup = True
quality = '3'
tmp_dir = 'tmp'

prefix = 'make_movie.py:  '
err = sys.stderr

#-------------------------------------------------------------------------------

def copyFiles(files, dir):
    """
        rename files to 'image????.EXT' where '????' are consecutive numbers
    """
    cnt = 0
    res = []
    for f in files:
        [main, ext] = os.path.splitext(f)
        name = dir + '/image%04i' % cnt + ext
        if not f == name:
            shutil.copyfile(f, name)
            #err.write("    %s --> %s\n" % ( f, name ))
        res.append(name)
        cnt += 1
    return res

#-------------------------------------------------------------------------------

def makeImages(format):
    """
        call `tool` to generate images in temporary directory
    """
    val = subprocess.call(tool + ['movie', 'image_dir='+tmp_dir, 'image_format='+format], stderr=None)
    if val:
        raise IOError("`%s' failed with value %i\n" % (tool[0], val))


def makeImagesUnzip(format, src='objects.cmo'):
    """
        unzip cytosim's input if necessary to generate the images
    """
    tmp_file = ''
    if not os.path.isfile(src):
        tmp_file = src + '.gz'
        if os.path.isfile(tmp_file):
            cmd = "gunzip %s -c > %s" % (tmp_file, src)
            subprocess.call(cmd, shell=True)
        else:
            tmp_file = ''
    if not os.path.isfile(src):
        raise IOError("file '%s' not found!\n" % src)
    makeImages(format)
    if tmp_file:
        os.remove(tmp_file)


def getImages(path, format):
    """
        Get images or call makeImage() to generate them
    """
    if tool and os.access(tool[0], os.X_OK):
        makeImagesUnzip(format)
        images = [os.path.join(tmp_dir, s) for s in os.listdir(tmp_dir)]
        images = sorted(images)
    else:
        # search for images already made:
        import glob
        images = sorted(glob.glob('m*[0-9].'+format))
        if not images:
            images = sorted(glob.glob('i*[0-9].'+format))
        err.write(prefix+"using %i images found inside `%s'\n" % (len(images), path))
    if not images:
        raise IOError("could not produce images!")
    return copyFiles(images, tmp_dir)


def getImageSize(file):
    """
        Call ffprobe, and parse output to extract size of video
    """
    res = [256, 256];
    proc = subprocess.Popen(['ffprobe', '-v', 'quiet', '-show_streams', file], stdout=subprocess.PIPE)
    if not proc.wait():
        for line in proc.stdout:
            [key, equal, value] = line.partition('=')
            if key=="width":
                res[0] = int(value)
            elif key=="height":
                res[1] = int(value)
    return res

#-------------------------------------------------------------------------------

def makeMovie(path, filename):
    """
        create movie.mp4 from PNG or PPM files in the current directory
        This entirely relies on ffmpeg
    """
    pattern = tmp_dir+'/image%04d.png'
    images = getImages(path, 'png')
    if not images:
        #print('looking for PPM images')
        pattern = tmp_dir+'/image%04d.ppm'
        images = getImages(path, 'ppm')
    if not images:
        raise IOError("no suitable images found")
    # build arguments for 'ffmpeg':
    args = ['ffmpeg', '-v', 'quiet', '-r', '%i'%rate, '-i', pattern]
    if codec == 'h265':
        args.extend(['-vcodec', 'libx265', '-tag:v', 'hvc1'])
    elif codec == 'h264':
        args.extend(['-vcodec', 'libx264'])
    else:
        args.extend(['-vcodec', 'mpeg4'])
    args.extend(['-pix_fmt', 'yuv420p', '-q:v', quality, filename])
    print(prefix + ' '.join(args))
    val = subprocess.call(args) #, stdout=None,stderr=None)
    if val:
        raise IOError("`ffmpeg` failed with value %i\n  %s\n" % (val, ' '.join(args)))

#-------------------------------------------------------------------------------

def process(path, filenames):
    """
        assemble movie in given directory (`path`)
    """
    global tmp_dir, cleanup
    res = 'movie.mp4'
    if ( res in filenames ) and lazy:
        err.write(prefix+"bailing out: `%s/%s' already exists!\n" % (path, res))
        return ''
    err.write("\n")
    try:
        import tempfile
        tmp_dir = tempfile.mkdtemp('', 'imgs-', '.')
        #err.write(prefix+"created directory %s\n" % tmp_dir)
    except Exception as e:
        err.write(prefix+"could not make temporary directory %s\n" % repr(e))
    try:
        makeMovie(path, res)
    except Exception as e:
        err.write(prefix+" %s\n" % str(e))
        return ''
    if cleanup:
        shutil.rmtree(tmp_dir)
        #err.write(prefix+"deleted directory %s\n" % tmp_dir)
    else:
        err.write(prefix+"folder `%s' contains generated images\n" % tmp_dir)
    return res


def process_dir(path):
    """call process() with appropriate arguments"""
    files = []
    with os.scandir() as it:
        for e in it:
            if e.is_file():
                files.append(e.name)
    return process(path, files)
    

#-------------------------------------------------------------------------------

def main(args):
    """
        process command line arguments
    """
    global tool, cleanup, lazy, codec, rate, quality
    paths = []
    
    try:
        arg = args[0]
    except:
        err.write(prefix+"you must specify an executable or a directory containing images\n")
        sys.exit()

    arg0 = os.path.expanduser(arg.split()[0])
    if os.access(arg0, os.X_OK) and os.path.isfile(arg0):
        tool = arg.split()
        tool[0] = os.path.abspath(arg0)
    elif os.path.isdir(arg):
        paths.append(arg)
    else:
        err.write(prefix+"You must specify an executable or a directory containing images\n")
        sys.exit()

    for arg in args[1:]:
        [key, equal, value] = arg.partition('=')
        
        if key=='' or equal!='=' or value=='':
            if os.path.isdir(arg):
                paths.append(arg)
            elif arg=='+':
                codec = 'h264'
            elif arg=='++':
                codec = 'h265'
            else:
                err.write(prefix+"ignored '%s' on command line\n" % arg)
        else:
            if key=='rate':
                rate = int(value)
            elif key=='cleanup':
                cleanup = bool(value)
            elif key=='codec':
                codec = value
            elif key=='quality':
                quality = value
            elif key=='lazy':
                lazy = int(value)
            else:
                err.write(prefix+"ignored '%s' on command line\n" % arg)

    if not paths:
        paths.append('.')

    cwd = os.getcwd()
    res = ''
    for p in paths:
        os.chdir(p)
        res = process_dir(p)
        os.chdir(cwd)
        # move file to parent directory if only one path was specified
        if res and len(paths) == 1 and not os.path.isfile(res):
            os.rename(os.path.join(paths[0],res), os.path.basename(res))
            err.write(prefix+"created %s\n" % res)
        else:
            err.write(prefix+"created %s/%s\n" % (p, res))

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])


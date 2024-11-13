// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include "frame_reader.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "filepath.h"
#include "messages.h"
#include "simul.h"
#include "dim.h"

// Use second definition to trace execution
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;

//------------------------------------------------------------------------------

FrameReader::FrameReader() : inputter(DIM)
{
    frameIndex = 0;
    lastLoaded = ~0UL;
}


void FrameReader::clear()
{
    inputter.clear();
}


void FrameReader::reset()
{
    inputter.rewind();
    clearPositions();
}


void FrameReader::openFile(std::string const& file)
{
    if ( !FilePath::is_file(file) && FilePath::is_file(file+".gz") )
    {
        std::string cmd = "gunzip " + file + ".gz";
        std::clog << cmd << '\n';
        if ( system(cmd.c_str()) )
            throw InvalidIO("failed to unzip `"+file+".gz'");
    }

    if ( inputter.open(file.c_str(), "rb") )
        return;
    
    if ( !inputter.file() )
        throw InvalidIO("file `"+file+"' not found");

    if ( inputter.error() )
        throw InvalidIO("file `"+file+"' is invalid");
 
    clearPositions();
    //std::clog << "FrameReader: has opened " << obj_file << '\n';
}


int FrameReader::badFile()
{
    if ( !inputter.file() )
        return 8;
    
    if ( inputter.eof() )
        inputter.clear();
    
    if ( !inputter.good() )
        return 7;
    
    return 0;
}

//------------------------------------------------------------------------------
#pragma mark -

void FrameReader::clearPositions()
{
    VLOG("FrameReader: clear\n");
    
    frameIndex = 0;
    framePos.clear();
    framePos.reserve(1024);
    framePos.resize(1);
    // store info for frame 0:
    fpos_t pos;
    if ( 0 == inputter.get_pos(pos) )
    {
        framePos[0].validity_ = 1;
        framePos[0].position_ = pos;
    }
}


void FrameReader::savePos(size_t frm, const fpos_t& pos, int confidence)
{
    if ( frm >= framePos.capacity() )
    {
        constexpr size_t chunk = 1024;
        size_t sz = ( frm + chunk - 1 ) & ~( chunk -1 );
        framePos.reserve(sz);
    }
    
    if ( frm >= framePos.size() )
    {
        size_t i = framePos.size();
        framePos.resize(frm+1);
        while ( i <= frm )
            framePos[i++].validity_ = 0;
    }
    
    if ( framePos[frm].validity_ < confidence )
    {
        framePos[frm].validity_ = confidence;
        framePos[frm].position_ = pos;
    
        //VLOG("FrameReader: position of frame " << frm << " is " << pos << '\n');
        VLOG("FrameReader: found position of frame "<<frm<<" ("<<confidence<<")\n");
    }
}


/**
 This uses the info stored in `framePos[]` to move to a position
 in the file where frame `frm` should start.
*/
size_t FrameReader::seekPos(size_t frm)
{
    if ( inputter.eof() )
        inputter.clear();
    
    if ( frm < 1 || framePos.empty() )
    {
        VLOG("FrameReader: seekPos rewind\n");
        inputter.rewind();
        return 0;
    }
    
    size_t inx = std::min(frm, framePos.size()-1);

    while ( inx > 0  &&  framePos[inx].validity_ == 0 )
        --inx;
    
    //check if we know already were the frame starts:
    if ( 0 < inx )
    {
        VLOG("FrameReader: using known position of frame " << inx << '\n');
        inputter.seek(framePos[inx].position_);
        return inx;
    }
    else {
        VLOG("FrameReader: rewind\n");
        inputter.rewind();
        return 0;
    }
}


size_t FrameReader::lastKnownFrame() const
{
    if ( framePos.empty() )
        return 0;
    size_t res = framePos.size()-1;
    while ( 0 < res  &&  framePos[res].validity_ < 2 )
        --res;
    return res;
}

//------------------------------------------------------------------------------
#pragma mark -


/// return code
enum FrameReaderCode { SUCCESS = 0, END_OF_FILE = 1, NOT_FOUND = 2, BAD_FILE = 4 };

/**
 scan file forward from current position to find the next Cytosim frame
 @return 0 if no frame was found
*/
int FrameReader::seekFrame(size_t frm)
{
    VLOG("FrameReader: seekFrame("<< frm <<")\n");
    
    size_t inx = seekPos(frm);
    
    if ( inx == frm )
        return SUCCESS;
    
    while ( ! inputter.eof() )
    {
        fpos_t pos;
        bool has_pos = false;
        size_t len = 1024;
        ssize_t read = 0;
        char * line = (char*)malloc(len);

        do {
            has_pos = !inputter.get_pos(pos);
            read = getline(&line, &len, inputter.file());

            if ( inputter.eof() )
                return END_OF_FILE;
            
#if 1 // backward compatibility code with format 42 before 2012
            if ( 7 < read && 0 == memcmp(line, "#frame ", 7) )
                break;
#endif
        } while ( read < 9 || memcmp(line, "#Cytosim ", 9) );
        
        //std::clog << "******\n";
        VLOG("           : " << line << '\n');
        free(line);

        if ( ! inputter.eof() )
        {
            if ( has_pos ) savePos(inx, pos, 2);
            if ( inx == frm )
            {
                if ( has_pos ) inputter.seek(pos);
                return SUCCESS;
            }
            ++inx;
        }
    }
    
    VLOG("FrameReader: seekFrame("<< frm <<") reached EOF\n");
    return END_OF_FILE;
}

//------------------------------------------------------------------------------
/**
 returns 0 for success, an error code, or throws an exception
 */
int FrameReader::loadFrame(Simul& sim, size_t frm, const bool reload)
{
    if ( badFile() )
        return BAD_FILE;

    VLOG("FrameReader: loadFrame("<<frm<<", reload="<<reload<<")\n");
    
    // what we are looking for might already be in the buffer:
    if ( frm == lastLoaded && ! reload )
        return SUCCESS;
    
    // it might be the next one in the buffer:
    if ( frm > 0  &&  frm-1 == lastLoaded )
        return loadNextFrame(sim);
    
    // otherwise, try to find the start tag from there:
    if ( SUCCESS != seekFrame(frm) )
        return NOT_FOUND;
    
    // store the position in the file:
    fpos_t pos;
    bool has_pos = !inputter.get_pos(pos);
    
    VLOG("FrameReader: reading file from "<<pos<<'\n');
    
    // read frame from file:
    int res = sim.reloadObjects(inputter);
    if ( !res )
    {
        VLOG("FrameReader: loadFrame("<<frm<<") successful\n");
        frameIndex = frm;
        lastLoaded = frameIndex;
        if ( has_pos )
            savePos(frameIndex, pos, 4);
        // the next frame should start at the current position:
        if ( 0 == inputter.get_pos(pos) )
            savePos(frameIndex+1, pos, 1);
        return SUCCESS;
    }
    else
    {
        VLOG("FrameReader: loadFrame("<<frm<<") failed: " <<res<< '\n');
        return res;
    }
}


/**
 returns 0 for success, an error code, or throws an exception
 */
int FrameReader::loadNextFrame(Simul& sim)
{
    if ( badFile() )
        return BAD_FILE;
    
    fpos_t pos;
    bool has_pos = !inputter.get_pos(pos);

    int res = sim.reloadObjects(inputter);
    if ( res == 0 )
    {
        if ( lastLoaded == frameIndex )
            ++frameIndex;
        lastLoaded = frameIndex;

        VLOG("FrameReader: loadNextFrame() read frame "<<currentFrame()<<'\n');
        
        // the position we used was good, to read this frame
        if ( has_pos )
            savePos(frameIndex, pos, 4);
        
        // the next frame should normally start from the current position:
        if ( !inputter.get_pos(pos) )
            savePos(frameIndex+1, pos, 1);
        return SUCCESS;
    }
    else
    {
        VLOG("FrameReader: loadNextFrame() failed at frame "<<currentFrame()<<'\n');
        return res;
    }
}


/**
 returns 0 for success, an error code, or throws an exception
 */
int FrameReader::loadLastFrame(Simul& sim, size_t cnt)
{
    if ( badFile() )
        return BAD_FILE;
    
    /// seek last known position:
    size_t frm = lastKnownFrame();
    if ( frm > 1 )
        inputter.seek(framePos[frm].position_);
    else
        inputter.rewind();
    
    /// go from here to last frame:
    int res = NOT_FOUND;
    while ( !sim.reloadObjects(inputter) )
    {
        frameIndex = frm++;
        lastLoaded = frameIndex;
        res = SUCCESS;
    }
    
    if ( res == SUCCESS && cnt > 0 )
    {
        frm = frm - 1 - cnt;
        // go back up by 'cnt' frames:
        if ( SUCCESS != seekFrame(frm) )
            return NOT_FOUND;
        
        if ( !sim.reloadObjects(inputter) )
            return NOT_FOUND;

        frameIndex = frm;
        lastLoaded = frameIndex;
        VLOG("FrameReader: loadFrame("<<frm<<") successful\n");
    }
    
    return res;
}

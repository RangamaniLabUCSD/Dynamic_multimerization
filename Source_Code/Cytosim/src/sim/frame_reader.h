// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "iowrapper.h"
#include <vector>

class Simul;


/// Helper class to access a particular frame in a trajectory file
/** 
 FrameReader is used to find a particular frame (eg. frame 10) in a trajectory file.
 FrameReader will handle basic IO failures and will remember the starting positions
 of all frames that were sucessfully identified, using this information to speed up
 future access to these and other frames in the same file. For example, if the 
 position of frame 50 is known and frame 60 is requested, it will start searching
 the file from the end of frame 50.

 FrameReader makes minimal assumptions on what constitutes a 'frame':
 - seekFrame() looks for a string that identifies the beggining of a frame (Cytosim).
 - loadFrame() calls Simul::reloadObjects() to read the content of the frame.
 .
 
 Frames are indexed starting at zero.
*/
class FrameReader
{
private:
    
    /// accessory class to store a position in a file
    class file_pos 
    {
    public:
        fpos_t position_; ///< starting position in the file
        int validity_;    ///< indicates if `position` is valid
        file_pos() { validity_ = 0; }
    };
    
    /// type for list of positions
    typedef std::vector<file_pos> PosList;
    
    /// the stream from which input is made
    Inputter inputter;
    
    /// starting position for each frame
    PosList framePos;
    
    /// index of frame stored currently
    size_t frameIndex;
    
    /// last frame loaded successfully
    size_t lastLoaded;
    
    /// go to a position where a frame close to `frm` is known to start
    size_t seekPos(size_t frm);

    /// remember position `pos` as the place where frame `frm` should start
    void savePos(size_t frm, const fpos_t& pos, int status);

public:
    
    /// constructor, after which openFile() should be called
    FrameReader();
    
    /// open file for input
    void openFile(std::string const& file);
    
    /// return index of current frame (0 is lowest valid value)
    size_t currentFrame() const { return frameIndex; }

    /// last frame seen in the file
    size_t lastKnownFrame() const;

    /// true when end of file is reached
    bool eof() const { return inputter.eof();  }
    
    /// true when end of file is reached
    bool good() const { return inputter.good();  }
    
    /// return 0 if file is good for input
    int badFile();

    /// rewind file
    void rewind() { inputter.rewind(); frameIndex=0; }
    
    /// dimensionality of vectors
    unsigned vectorSize() const { return inputter.vectorSize(); }
    
    /// clear error flag of file
    void clear();
    
    /// clear the record of positions of frames in file
    void clearPositions();

    /// rewind file and clear position buffer
    void reset();
    
    /// find the starting point of frame `frm` by brute force !
    int seekFrame(size_t frm);
    
    /// load specified frame into given Simul (the index of first frame is 1)
    int loadFrame(Simul&, size_t frm, bool reload = false);
    
    /// read the next frame in the file, return 0 for SUCCESS, 1 for EOF
    int loadNextFrame(Simul&);
    
    /// read the last frame in the file, return 0 for SUCCESS, 1 if no frame was found
    int loadLastFrame(Simul&, size_t cnt = 0);

};



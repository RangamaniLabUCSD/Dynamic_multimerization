// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "sim.h"
#include "meca.h"
#include "modulo.h"
#include "glossary.h"
#include "simul_prop.h"
#include "print_color.h"
#include "display1.h"
#include "display2.h"
#include "display3.h"
#include "save_image.h"
#include "filepath.h"
#include <unistd.h>
#include <cstdlib>
#include "glut.h"

extern void helpKeys(std::ostream&);

extern Modulo const* modulo;

//------------------------------------------------------------------------------
#pragma mark -

void Player::setStyle(const unsigned style)
{
    if ( mDisplay )
    {
        //restore the previous OpenGL state
        glPopAttrib();
        delete(mDisplay);
        mDisplay = nullptr;
    }
    
    //save the current OpenGL state
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    
    switch ( style )
    {
        default:
        case 1: mDisplay = new Display1(&disp);  break;
        case 2: mDisplay = new Display2(&disp);  break;
        case 3: mDisplay = new Display3(&disp);  break;
    }
    disp.style = style;

    //initialize Views associated with opened GLUT windows:
    for ( size_t n = 1; n < glApp::views.size(); ++n )
    {
        View & view = glApp::views[n];
        if ( view.window() > 0 )
        {
            //std::clog << "initializing GLUT window " << n << std::endl;
            //glutSetWindow(view.window());
            view.initGL();
            glViewport(0, 0, view.width(), view.height());
        }
    }
}


/**
 Build a message containing the label and the time.
 In live mode, it also adds 'Live' or it indicates the frame index,
 and the force generated by the mouse-controlled Single.
 */
std::string Player::buildLabel() const
{
    std::ostringstream oss;
    oss.precision(3);

    oss << std::setw(8) << std::fixed << simul.time() << "s";
    
    //display the force exerted by the mouse-controled Single:
    Single const* sh = thread.handle();
    if ( sh && sh->attached() )
        oss << "\nHandle: " << sh->force().norm() << "pN";

    if ( thread.alive() && goLive )
    {
        oss << "\nLive";
        //display ratio number-of-time-step / frame
        if ( prop.period > 1 )
            oss << " " << prop.period;
    }
    else if ( thread.currentFrame() > 0 )
    {
        oss << "\nFrame " << thread.currentFrame();
    }

    return oss.str();
}


/**
 This information is displayed in the top corner of the window.
 Calling simul.report() maked sure that the message is identical to what
 would be printed by the command 'report'.
 */
std::string Player::buildReport(std::string arg) const
{
    if ( ! arg.empty() )
    {
        Glossary glos;
        // separate options:
        std::string::size_type pos = arg.find(' ');
        if ( pos != std::string::npos )
        {
            glos.read_string(arg.substr(pos+1).c_str(), 2);
            arg = arg.substr(0, pos);
        }
        try
        {
            std::stringstream ss;
            simul.report(ss, arg, glos);
            std::string res = ss.str();
            if ( res.size() > 1  &&  res.at(0) == '\n' )
                return res.substr(1);
            return res;
        }
        catch ( Exception & e )
        {
            return e.what();
        }
    }
    return "";
}

/**
 This text is normally displayed in the center of the window
 */
std::string Player::buildMemo(int type) const
{
    std::ostringstream oss;
    switch ( type )
    {
        case 0: return "";
        case 1: return "Please, visit www.cytosim.org";
        case 2: helpKeys(oss); return oss.str();
        case 3: glApp::help(oss);  return oss.str();
        case 4: writePlayParameters(oss, true); return oss.str();
        case 5: writeDisplayParameters(oss, true); return oss.str();
    }
    return "";
}

//------------------------------------------------------------------------------
#pragma mark - Display

void Player::autoTrack(FiberSet const& fibers, View& view)
{
    Vector G(0, 0, 0);
    Vector D(0, 0, 0);
    real vec[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    
    if ( view.track_fibers & 1 )
    {
        Vector M, P;
        FiberSet::infoPosition(fibers.collect(), M, G, P);
        view.move_shift(Vector3(G));
        //std::clog << "auto center: " << G << std::endl;
    }
    
    if ( view.track_fibers & 2 )
    {
        // align with mean nematic direction
        FiberSet::infoNematic(fibers.collect(), vec);
        view.align_with(Vector3(vec));
        //view.rotation.setFromMatrix3(vec);
        //view.rotation.conjugate();
        //std::clog << "auto rotate: " << Vector3(vec) << std::endl;
    }

    if ( view.track_fibers & 4 )
    {
        real sum = 0;
        real avg[3] = { 0 };
        real mom[9] = { 0 };
        FiberSet::infoComponents(fibers.collect(), sum, avg, mom, vec);
        // get rotation from matrix:
        view.rotation.setFromMatrix3(vec);
        // inverse rotation:
        view.rotation.conjugate();
        //std::clog << "auto quat: " << view.rotation << std::endl;
    }
}


/**
 Adjust to see the biggest Space in simul
 */
void Player::autoScale(SpaceSet const& spaces, View& view)
{
    real rad = 0;
    for ( Space const* spc = spaces.first(); spc; spc=spc->next() )
        rad = std::max(rad, spc->max_extension());
    if ( rad > 0 )
    {
        //std::clog << "auto_scale " << rad << '\n';
        view.view_size = GLfloat(2*rad);
        view.zoom_in(0.933033);
        --view.auto_scale;
    }
}


void Player::prepareDisplay(View& view, int mag)
{    
    //gle::gleReportErrors(stderr, "before prepareDisplay");
    
    //----------------- automatic adjustment of viewing area:

    if ( view.auto_scale > 0 )
        autoScale(simul.spaces, view);
    
    //----------------- auto-track:
    
    if ( view.track_fibers )
        autoTrack(simul.fibers, view);
    
    //----------------- texts:
    
    view.setLabel(buildLabel());
    view.setMessage(buildReport(prop.report));
    
    //----------------- set pixel size and unit-size:
    /*
     if disp.point_value is set, line-width and point-size were specified in 'real' units,
     and otherwise, they were specified in pixels.
     */

    GLfloat pix = view.pixelSize();
    //std::clog << " pixel size = " << pix << '\n';

    if ( disp.point_value > 0 )
        mDisplay->setPixelFactors(pix/mag, mag*disp.point_value/pix);
    else
        mDisplay->setPixelFactors(pix/mag, mag);

    gle::gleReportErrors(stderr, "before prepareDisplay");

    try {
        mDisplay->setStencil(view.stencil && ( DIM == 3 ));
        mDisplay->prepareForDisplay(simul, dproperties);
        //std::clog << " dproperties.size() = " << dproperties.size() << '\n';
    }
    catch(Exception & e) {
        std::cerr<<"Error in prepareDisplay: " << e.what() << '\n';
    }
}

//------------------------------------------------------------------------------
void Player::displayCytosim()
{
    // clear pixels:
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    try {
        // draw:
        if ( modulo && disp.tile )
            mDisplay->displayTiled(simul, disp.tile);
        else
            mDisplay->display(simul);

        if ( disp.draw_links )
        {
#if DRAW_MECA_LINKS
            glPushAttrib(GL_LIGHTING_BIT);
            glDisable(GL_LIGHTING);
            glLineWidth(4);
            glPointSize(8);
            glEnable(GL_LINE_STIPPLE);
            simul.drawLinks();
            glDisable(GL_LINE_STIPPLE);
            glPopAttrib();
            gle::gleReportErrors(stderr, "Simul::drawLinks()");
#endif
        }
    }
    catch(Exception & e) {
        std::cerr << "Error in display: " << e.what() << std::endl;
    }
}


void Player::readDisplayString(View& view, std::string const& str)
{
    //std::clog << "readDisplayString " << str << std::endl;
    try
    {
        Glossary glos(str);
        disp.read(glos);
        const int W = view.width(), H = view.height();
        view.read(glos);
        // window size cannot be changed:
        view.window_size[0] = W;
        view.window_size[1] = H;
    }
    catch( Exception & e )
    {
        std::cerr << "Error while reading simul:display: " << e.what();
    }
}


/**
 This display the full Scene
 */
void Player::displayScene(View& view, int mag)
{
    if ( simul.prop->display_fresh )
    {
        readDisplayString(view, simul.prop->display);
        simul.prop->display_fresh = false;
    }
    //thread.debug("display");
    prepareDisplay(view, mag);
    view.openDisplay();
    displayCytosim();
    view.closeDisplay();
    glFinish();
}

//------------------------------------------------------------------------------
#pragma mark - Export Image

/**
 Export image from the current OpenGL back buffer,
 in the format specified by 'PlayerProp::image_format',
 in the current working directory
 */
int Player::saveView0(const char* filename, const char* format, int downsample) const
{
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);
    int err = SaveImage::saveImage(filename, format, vp, downsample);
    if ( err == 0 && simul.prop->verbose > 0 )
    {
        printf("\r saved %ix%i (%i) snapshot %s    ", vp[2]/downsample, vp[3]/downsample, downsample, filename);
        fflush(stdout);
    }
    return err;
}


/**
 Export image from the current OpenGL back buffer,
 in the format specified by 'PlayerProp::image_format',
 in the folder specified in `PlayerProp::image_dir`.
 The name of the file is formed by concatenating 'root' and 'indx'.
 */
int Player::saveView(const char* root, size_t indx, int downsample) const
{
    char const* fmt = prop.image_format.c_str();
    char str[1024] = { 0 };
    snprintf(str, sizeof(str), "%s%04lu.%s", root, indx, fmt);
    int cwd = FilePath::change_dir(prop.image_dir, true);
    int err = saveView0(str, fmt, downsample);
    FilePath::change_dir(cwd);
    return err;
}

//------------------------------------------------------------------------------

void displayMagnified(int mag, void * arg)
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    static_cast<Player*>(arg)->displayCytosim();
}


/**
 save an image where the resolution is magnified by a factor `mag`.
 This requires access to the simulation world.
 */
int Player::saveScene(const int mag, const char* name, const char* format, const int downsample)
{
    View & view = glApp::currentView();
    const int W = view.width(), H = view.height();
    thread.lock();
    
    //std::clog << "saveMagnifiedImage " << W << "x" << H << " mag=" << mag << '\n';
    prepareDisplay(view, mag);
    view.openDisplay();
    int err = SaveImage::saveMagnifiedImage(mag, name, format, W, H, displayMagnified, this, downsample);
    if ( err )
        err = SaveImage::saveCompositeImage(mag, name, format, W, H, view.pixelSize(), displayMagnified, this, downsample);
    if ( !err )
        printf("saved %ix%i snapshot %s\n", mag*W/downsample, mag*H/downsample, name);
    view.closeDisplay();
    
    thread.unlock();
    return err;
}


/**
 save an image where the resolution is magnified by a factor `mag`.
 This requires access to the simulation world.
 */
int Player::saveScene(const int mag, const char* root, unsigned indx, const int downsample)
{
    char str[1024] = { 0 };
    char const* format = prop.image_format.c_str();
    snprintf(str, sizeof(str), "%s%04i.%s", root, indx, format);
    int cwd = FilePath::change_dir(prop.image_dir, true);
    int err = saveScene(mag, str, format, downsample);
    FilePath::change_dir(cwd);
    glApp::postRedisplay();
    return err;
}

// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dimerizer.h"
#include "dimerizer_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"


Dimerizer::Dimerizer(DimerizerProp const* p, HandMonitor* h)
: Hand(p,h), prop(p)
{
}

//detach
void Dimerizer::detach()
{
    assert_true (!fbFiber)
    auto ha2=haHand;
    // std::cout<<"Hands detached "<<this<<" "<<ha2<<std::endl;
    // std::cout<<"nextDetach "<<nextDetach<<std::endl;
    assert_true( attached() );
    haMonitor->beforeDetachment(this);
    haHand = nullptr;
    marktodetach = false;
    hasbeentossed = false;
}

/**
tests detachment
 */
void Dimerizer::stepUnloaded()
{
    assert_true( attached() );
    testDetachment();
}

void Dimerizer::stepLoaded(Vector const& force, real force_norm)
{
    //Only do unbind coin tosses from one of the two hands.
    // if (this<haHand)
    // {
    assert_true( attached() );
    //Only one of the two hands involved in the dimerizer have
    //nextDetach. So we need to find use that. That will be used 
    //in the testKramersDetachment function;
    //assert_true( nextDetach >= 0 );

    testKramersDetachment(force_norm);
    // }
}

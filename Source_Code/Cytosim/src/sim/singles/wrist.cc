// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "wrist.h"
#include "simul.h"
#include "meca.h"
#include "modulo.h"
#include "dimerizer.h"


extern Modulo const* modulo;


Wrist::Wrist(SingleProp const* sp, Mecable const* mec, const unsigned pti)
: Single(sp)
{
    anchor.set(mec, pti);

#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter("single:diffusion cannot be > 0 if activity=anchored");
#endif
}


Wrist::Wrist(SingleProp const* sp, Mecable const* mec, unsigned a, unsigned b, real c)
: Single(sp)
{
    assert_true(mec);
    anchor.set(mec, a, b, c);
    
#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter("single:diffusion cannot be > 0 if activity=anchored");
#endif
}


Wrist::Wrist(SingleProp const* sp, Mecable const* mec, unsigned ref, Vector pos)
: Single(sp)
{
    assert_true(mec);
    anchor.set(mec, ref, pos);

    
#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter("single:diffusion cannot be > 0 if activity=anchored");
#endif
}


Wrist::~Wrist()
{
}


Vector Wrist::force() const
{
    assert_true( sHand->attached() );
    Vector d = posFoot() - sHand->pos();
    
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


bool Wrist::stepF(Simul& sim)
{
    // std::cout<<"stepF Hand "<<sHand<<" Hand "<<sHand->hand()<<" Fiber "<<sHand->fiber()<<" TAG "<<sHand->tag()<<std::endl;
    // if(sHand->tag()==Dimerizer::TAG){
    //     Dimerizer* d2 = static_cast<Dimerizer*>(sHand);
    // std::cout<<"hasbeentossed "<<d2->hasbeentossed<<" marktodetach "<<d2->marktodetach<<std::endl;
    // }

    //if the hand is bound and sHand is NOT nullptr, return
    // if(sHand->tag()==Dimerizer::TAG &&
    // sHand->attached() == true && sHand->hand()!=nullptr)
    //     return;
    assert_false( sHand->attached() );
    //std::cout<<"Hand attached status "<<sHand->attached() <<std::endl;
    return sHand->stepUnattached(sim, posFoot());
}


void Wrist::stepA()
{
    // std::cout<<"stepA Hand "<<sHand<<" attached hand "<<sHand->hand()<<" fiber "<<sHand->fiber()<<std::endl;
    // std::cout<<"attached? "<<sHand->attached()<<" null? "<<(sHand->hand()==nullptr)<<" TAG check "<<
    // (sHand->tag())<<" binding key "<<sHand->prop->binding_key<<std::endl;
    //if the hand is free and sHand is nullptr, return
    // if(sHand->tag()==Dimerizer::TAG && sHand->attached() == false && sHand->hand()==nullptr)
    //     return;
    
    assert_true( sHand->attached() );
    Vector f = force();
    sHand->stepLoaded(f, f.norm());
    // std::cout<<"attached? "<<sHand->attached()<<" null? "<<(sHand->hand()==nullptr)<<" TAG check "<<
    // (sHand->tag())<<" binding key "<<sHand->prop->binding_key<<std::endl;
}

void    Wrist::steplastF()
{
    ///check if hand is dimerizer
    // std::cout<<"last Hand "<<sHand<<std::endl;
    if (sHand->tag()==Dimerizer::TAG){
        auto otherhand = sHand->hand();
        //std::cout<<"tmphaHand exists? "<<otherhand->tmphaHand<<std::endl;
        ///check if tmphaHand exists
        if( otherhand->tmphaHand ){
            assert_true(otherhand->tmphaHand==sHand);
            assert_false(otherhand->attached());
            otherhand->attach(sHand);
            otherhand->tmphaHand=nullptr;
        }
        // std::cout<<"tmphaHand exists? "<<sHand->tmphaHand<<std::endl;
        // ///check if tmphaHand exists
        // if( sHand->tmphaHand ){
        //     assert_false(attached())
        //     sHand->attach(sHand->tmphaHand);
        //     sHand->tmphaHand=nullptr;
        // }
    }
}


void Wrist::setInteractions(Meca & meca) const
{
    //Need to edit to add interaction here.
    //See chain.cc line 1421
    if (sHand->tag() == Dimerizer::TAG){
        ////////////////--thishand->boundhand->Single(Wrist)
        Wrist* wrist2 = static_cast<Wrist*>(sHand->hand()->handmonitor());
        anchor.addLink(meca, wrist2->mecapoint(), prop->stiffness);
    }
    else
        anchor.addLink(meca, sHand->interpolation(), prop->stiffness);
}


void Wrist::write(Outputter& out) const
{
    sHand->write(out);
    out.writeSoftSpace();
    anchor.write(out);
}


void Wrist::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    const int s = state();

    sHand->read(in, sim);
    
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 47 )
    {
        Mecapoint base;
        base.read(in, sim);
        anchor.set(base.mecable(), base.point());
    }
    else
#endif
        anchor.read(in, sim);
    
    /*
     Because the SingleSet has 2 lists where Single are stored depending
     on their bound/unbound state, we need to unlink and relink here, in
     case the state stored on file is different from the current state.
     */
    if ( s != state() )
    {
        SingleSet * set = static_cast<SingleSet*>(objset());
        if ( set )
        {
            if ( s )
                set->relinkD(this);
            else
                set->relinkA(this);
        }
    }
}


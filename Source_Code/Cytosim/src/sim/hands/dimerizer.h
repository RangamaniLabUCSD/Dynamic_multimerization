// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DIMERIZER_H
#define DIMERIZER_H

#include "hand.h"
#include "dimerizer_prop.h"

/// a Hand that can move smoothly on a Fiber
/**
 The Dimerizer is a Hand, and thus can bind and unbind from fibers.
 
 A bound Dimerizer can move along its fiber.
 The direction of the dimerizer is set by the sign of @ref DimerizerPar "unloaded_speed".
 If the speed is positive, the dimerizer attempts to move towards the PLUS_END.
 
 The speed is linearly proportional to the load of the dimerizer.
 The load is the projection of the force vector on the direction of the fiber.
 
     real load = force * direction_of_fiber;
     real speed = unloaded_speed * ( 1 + load / stall_force );
 
 The actual movement depends on the time step:

     real displacement = speed * time_step;
 
 As defined in Hand, detachment increases exponentially with force.

 See Examples and the @ref DimerizerPar.
 @ingroup HandGroup
 */
class Dimerizer : public Hand
{
private:
    
    /// disabled default constructor
    Dimerizer();

public:
    /// Object::TAG = 'v' represents the 'void' pointer
    static const HandTag TAG = 'd';
    /// a character identifying the class of this object
    virtual HandTag tag() const { return TAG; }

    /// Used to ensure poisson sampling for unbinding only happens once
    bool hasbeentossed = false;
    bool marktodetach = false;
    
    /// Property
    DimerizerProp const* prop;
    
    /// constructor
    Dimerizer(DimerizerProp const*, HandMonitor*);
    
    /// destructor
    ~Dimerizer() {}

    ///detach
    virtual void detach();
    /**
     Test for spontaneous detachment using Gillespie approach.
     @return true if the test has passed, and detach() was called.
     see @ref Stochastic
     */
    bool testDetachment()
    {
        nextDetach -= prop->unbinding_rate_dt;
        
        if ( nextDetach <= 0 )
        {
            if(haHand->tag()==Dimerizer::TAG) haHand->detach();
            detach();
            return true;
        }
        
        return false;
    }
    
    
    /**
     Test for spontaneous detachment using Gillespie approach.
     @return true if the test has passed, and detach() was called.
     see @ref Stochastic
     */
    // Here, we would have this specialization for dimerizer.
    bool testKramersDetachment(const real force)
    {       
        ///Make sure this hand has not tried to unbind
        assert_false(hasbeentossed);
        /// List of various permutations
        /// |d2->hasbeentossed|marktodetach|Comment
        /// |true*            |true        |this->detach();
        /// |false            |false       |set this as tossed, try to unbind
        /// |                 |            |success:detach, mark the other to detach
        /// |                 |            |failure: return
        /// |true             |false       |reset d2->hasbeentossed=false; return
        /// |false            |true        |logic does not allow**
        /// *Note: if the other hand detaches,reset hapens: d2->hasbeentossed= false 
        /// **Note: Even though logically this is not allowed, this is how detachment
        /// happens as reset hapens when other hand detaches: d2->hasbeentossed= false
        /// check if this hand has been marked to detach
        if(marktodetach){
            detach();
            return false;
        }
        Dimerizer* d2=static_cast<Dimerizer*>(haHand);
        // std::cout<<"Hand 1 hasbeentossed "<<hasbeentossed<<" marktodetach "<<marktodetach<<std::endl;
        // std::cout<<"Hand 2 hasbeentossed "<<d2->hasbeentossed<<" marktodetach "<<d2->marktodetach<<std::endl;
        real ndtemp;
        /// check if the other hand has tried to unbind
        if(d2->hasbeentossed) {
            //reset
            d2->hasbeentossed = false;
            return false;
        }
        else{
            hasbeentossed = true;
            real otherhand_nextDetach = d2->getnextDetach();
            //If this is the hand with the accurate nextDetach
            if ( nextDetach >=0 ){
                nextDetach -= prop->unbinding_rate_dt * exp(force*prop->unbinding_force_inv);
                ndtemp = nextDetach;
            }
            //If the other hand has that info
            else{
                ndtemp = otherhand_nextDetach;
                ndtemp -= prop->unbinding_rate_dt * exp(force*prop->unbinding_force_inv);
                d2->setnextDetach(ndtemp);
            }
            /*
            Attention: the exponential term can easily become numerically "infinite",
            which is problematic if 'unbinding_rate==0' and 'unbinding_force' is finite.
            This issue is handled in HandProp::complete()
            */
            if ( ndtemp <= 0 ){
                /// mark other hand to be detached
                d2->marktodetach = true;
                detach();
                return true;
            }
        }
    return false;
    }

    /// simulate when `this` is attached but not under load
    void   stepUnloaded();
    
    /// simulate when `this` is attached and under load
    void   stepLoaded(Vector const& force, real force_norm);
    
    ///getter and setter for nextDetach
    void setnextDetach(real n) {nextDetach = n;}
    real getnextDetach() {return nextDetach;}
};

#endif


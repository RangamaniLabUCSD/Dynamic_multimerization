// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DIMERIZER_PROP_H
#define DIMERIZER_PROP_H

#include "hand_prop.h"

/// Additional Property for Dimerizer
/**
 @ingroup Properties
 */
class DimerizerProp : public HandProp
{
    friend class Dimerizer;
    
public:
    
    /**
     @defgroup DimerizerPar Parameters of Dimerizer
     @ingroup Parameters
     Inherits @ref HandPar.
     @{
    */
    
    /// force at which speed reaches zero (must be > 0)
    /**
     For kinesin, see for example:
     http://dx.doi.org/10.1016/0092-8674(94)90060-4
     http://www.ncbi.nlm.nih.gov/pmc/articles/PMC42784/
     For dynein, see:
     
     For myosin, see:
     http://dx.doi.org/10.1038/368113a0
     */
    real    stall_force;
    
    
    /// speed of the dimerizer when its load is zero
    /** 
     A positive value specifies a plus-end directed dimerizer.
     A negative value specifies a minus-end directed dimerizer.
     */
    real    unloaded_speed;
    
    /// if true, the speed is limited to the range [0, 2*unloaded_speed]
    /** 
     With ( limit_speed = 1 ), a plus-end directed dimerizer will never move towards
     the minus-end, and will not exceed 2x its unloaded speed.
     For a minus-end directed dimerizer, the permitted range is [2*unloaded_speed, 0],
     thus also excluding backward steps and excessive speed.
     */
    bool    limit_speed;
    
    /// @}
    
private:
    
    /// limits for a displacement in one time_step apply if ( limit_speed = true )
    real    min_dab, max_dab;

    /// variables derived from `unloaded_speed`
    real    set_speed_dt, abs_speed_dt, var_speed_dt;
    
public:

    /// constructor
    DimerizerProp(const std::string& n) : HandProp(n)  { clear(); }
    
    /// destructor
    ~DimerizerProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// perform additional tests for the validity of parameters, given the elasticity
    void checkStiffness(real stiff, real len, real mul, real kT) const;
    
    /// return a carbon copy of object
    Property* clone() const { return new DimerizerProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif


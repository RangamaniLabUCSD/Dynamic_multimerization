// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "common.h"
#include "property_list.h"
#include "dimerizer_prop.h"
#include "dimerizer.h"
#include "simul.h"


Hand * DimerizerProp::newHand(HandMonitor* m) const
{
    return new Dimerizer(this, m);
}


void DimerizerProp::clear()
{
    HandProp::clear();
    
    stall_force       = 0;
    unloaded_speed    = 0;
    limit_speed       = true;
    var_speed_dt = 0;
    set_speed_dt = 0;
    abs_speed_dt = 0;
}


void DimerizerProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    glos.set(stall_force,    "stall_force")    || glos.set(stall_force,    "force");
    glos.set(unloaded_speed, "unloaded_speed") || glos.set(unloaded_speed, "speed");
#ifdef BACKWARD_COMPATIBILITY
    glos.set(unloaded_speed, "max_speed");
#endif
    glos.set(limit_speed,    "limit_speed");
}


void DimerizerProp::complete(Simul const& sim)
{
    HandProp::complete(sim);
    
    // if ( sim.ready() && stall_force <= 0 )
    //     throw InvalidParameter("dimerizer:stall_force must be > 0");
    
    // set_speed_dt = sim.time_step() * unloaded_speed;
    // abs_speed_dt = fabs(set_speed_dt);
    // var_speed_dt = abs_speed_dt / stall_force;
    
    // // The limits for a displacement in one time_step apply if ( limit_speed = true )
    // if ( unloaded_speed > 0 )
    // {
    //     min_dab = 0;
    //     max_dab = 2 * sim.time_step() * unloaded_speed;
    // }
    // else
    // {
    //     min_dab = 2 * sim.time_step() * unloaded_speed;
    //     max_dab = 0;
    // }
}


void DimerizerProp::checkStiffness(real stiff, real len, real mul, real kT) const
{
    HandProp::checkStiffness(stiff, len, mul, kT);
}


void DimerizerProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "stall_force",       stall_force);
    write_value(os, "unloaded_speed",    unloaded_speed);
    write_value(os, "limit_speed",       limit_speed);
}


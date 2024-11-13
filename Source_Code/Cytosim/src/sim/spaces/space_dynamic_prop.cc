// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cmath>
#include "space_dynamic_prop.h"
#include "space_dynamic_ellipse.h"
#include "space_prop.h"
#include "property_list.h"
#include "glossary.h"


//------------------------------------------------------------------------------


Space * SpaceDynamicProp::newSpace() const
{
	return new SpaceDynamicEllipse(this);
}


void SpaceDynamicProp::clear()
{
    tension = 0;
    volume  = 0;
	SpaceProp::clear();
}

void SpaceDynamicProp::read(Glossary& glos)
{
    SpaceProp::read(glos);
	glos.set(tension, "tension");
}



void SpaceDynamicProp::complete(Simul const& sim)
{
	SpaceProp::complete(sim);
	
	if	(tension < 0)
		throw InvalidParameter("tension must be positive");
}


void SpaceDynamicProp::write_values(std::ostream& os) const
{
    SpaceProp::write_values(os);
	write_value(os, "tension", tension);
}

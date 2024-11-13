# Tutorial: Chromosome Capture

This is unfinished!

Author: XXX (XX.09.2018)


### Objective

This tutorial will teach you how to simulate the capture of chromosomes by an aster of dynamic microtubules.

The configuration is directly taken from:  
>*[F-Actin nucleated on chromosomes coordinates their capture by microtubules in oocyte meiosis
](http://doi.org/10.1083/jcb.201802080)*  
>Mariia Burdyniuk, Andrea Callegari, Masashi Mori, François Nédélec, and Péter Lénárt.  
>Journal of Cell Biology 2018.

We will start with [this configuration file](data/capture.cym).
Let's see what this files contains.

# Code

Important: this simulation used two custom features, that are disabled by default in the code:

	NEW_LENGTH_DEPENDENT_CATASTROPHE

	NEW_RADIAL_FLOW

These compile switches should be set to 1, and the code needs to be recompiled: `make clean; make`

	#define NEW_LENGTH_DEPENDENT_CATASTROPHE 1

	#define NEW_RADIAL_FLOW 1


# Setup

	
	set simul system
	{
	    steric = 0, 500;
	    time_step = 0.05;
	    viscosity = 0.1;
	    display = ( style = 1; );
	}
	
	set space cell
	{
	    shape = sphere;
	}
		
	new cell
	{
	    radius = 35;
	}


# Microtubules

	set fiber microtubule
	{
	    activity = classic;
	    catastrophe_length = 37;
	    catastrophe_outside = 1;
	    catastrophe_rate = 0.15;
	    confine = inside, 100;
	    persistent = 1;
	    growing_force = 1.67;
	    growing_speed = 0.5;
	    min_length = 0.25;
	    rebirth_rate = 1;
	    rescue_rate = 0;
	    rigidity = 22;
	    segmentation = 3;
	    shrinking_speed = -1;
	}
	
	set hand nucleator
	{
	    activity = nucleate;
	    nucleate = 1, microtubule, ( length = 3; plus_end = grow; );
	    unbinding = 0, 3;
	    display = ( color = green; size = 5; );
	}
	
	set single complex
	{
	    activity = fixed;
	    hand = nucleator;
	    stiffness = 1000;
	}
	
	new 500 complex
	{
	    position = sphere 2 at 4 32 0;
	}
	
	new 500 complex
	{
	    position = sphere 2 at -4 32 0;
	}

# Kinetochores (capturing entities)
	
	set hand glue
	{
	    activity = move;
	    binding_range = 0.25;
	    binding_rate = 0;
	    hold_shrinking_end = 1;
	    unloaded_speed = -0.33;
	    stall_force = 1;
	    unbinding = 0, inf;
	    display = ( color = green; size = 9; );
	}
	
	set single kinetochore
	{
	    hand = glue;
	    stiffness = 100;
	}
	
# Chromosomes

	
	set solid chromosome
	{
	    confine = inside, 100;
	    flow_center = 0 33 0;
	    flow_time = 750, 900;
	    steric = 1;
	    display = ( color = 0xFFFFFF44; style = 7; );
	}

	new 22 chromosome
	{
	    point1 = center, 0.8;
	    point2 = +0.8 0 0, 0, kinetochore;
	    point3 = -0.8 0 0, 0, kinetochore;
	    position = sphere 35;
	}
	
# Simulation
	
	run system
	{
	    nb_frames = 120;
	    nb_steps = 4800;
	}
	
	%%%%%%%%%%%%%%%%%%%% STAGE 1
	
	change glue
	{
	    binding_rate = 0.1;
	}
	
	run system
	{
	    nb_frames = 30;
	    nb_steps = 1200;
	}
	
	%%%%%%%%%%%%%%%%%%%% STAGE 2
	
	change glue
	{
	    binding_rate = 0.3;
	}
	
	run system
	{
	    nb_frames = 30;
	    nb_steps = 1200;
	}
	
	%%%%%%%%%%%%%%%%%%%% STAGE 3
	
	change glue
	{
	    binding_rate = 0.7;
	}
	
	run system
	{
	    nb_frames = 30;
	    nb_steps = 1200;
	}
	
	%%%%%%%%%%%%%%%%%%%% STAGE 4
	
	change glue
	{
	    binding_rate = 1;
	}
		
	run system
	{
	    nb_frames = 240;
	    nb_steps = 9600;
	}
	

### The end!

Congratulation, you have completed the tutorial.

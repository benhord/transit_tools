import numpy as np
import matplotlib.pyplot as plt
import rebound

#function to perform a dynamical simulation with REBOUND
def plsim(mstar, plparams=None, plmass=None, a=None, eflag=False):
    #!!Add ability to gather stellar params on its own or just let main do
    #  that part?!!
    """
    Function to perform a dynamical simulation with the REBOUND code given 
    arbitrary planet inputs. Order of inputs for each planet parameters does not
    matter, as long as the order is consistent between each planet parameter.
    Ordering by distance from center of system is recommended.

    Parameters
    ----------
    mstar : float
        The mass of the central object in the simulation. Units in solar masses.
    plparams : list of dicts
        List of dictionaries of planet parameters. Each planet should have its
        own dictionary in the list. Mutually exclusive with specifying these
        parameters individually through other arguments.
    plmass : list or float
        List of planet masses in order of distance from center of system. 
        Alternatively, a float can be specified if there is only one planet.
        Mututally exclusive with plparams.
    a : list or float
        Semi-major axis of planets in order of distance from central body in 
    eflag : bool
        Flag determining whether eccentricity will be simulated. May make the
        simulation run slower.

    Returns
    -------
    """
    #set up simulation and add central star
    sim = rebound.Simulation()
    sim.add(m=mstar)

    #check inputs for planet parameters
    

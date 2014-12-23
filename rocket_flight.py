#!/usr/bin/env python2

"""
Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
Time-stamp: <2014-09-16 00:28:56 (jonah)>

This is my solution to the rocket flight problem.

Strategy:

We evolve a system of three variables, h, v, and mp, where:
h  = the altitude of the rocket
v  = the vertical velocity of the rocket
mp = the mass of the propellant remaining on the rocket.
"""

# Imports
# ------------------------------------------------------------------
import numpy
import matplotlib.pyplot as plt
# ------------------------------------------------------------------


# Set the matplotlib rcparams
# -----------------------------------------------------------------
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16
# -----------------------------------------------------------------


# Global constants
# ------------------------------------------------------------------
ms = 50 # kg            -- Weight of the rocket shell
g = 9.81 # m/s          -- Acceleration due to gravity
rho = 1.091 # kg/m^3    -- Average air density
r = 0.5 # m             -- Radius of rocket
A = numpy.pi*r**2 # m^2 -- Cross-sectional area of rocket
ve = 325 # m/s          -- Exhaust speed
CD = 0.15 #             -- Drag coefficient
mp0 = 100 # kg          -- Initial weight of the rocket propellant
dt = 0.1 # s            -- time step

# Burn rate of the propellant. Not really a constant, but
# conceptually the same.
burn_rate = lambda time: 20 if time < 5 else 0 # kg/s

# Initial data we choose. Assuming initial height and velocity are zero
# but this isn't really necessary.
h0 = 0.
v0 = 0.
# ------------------------------------------------------------------


# Functions
# ------------------------------------------------------------------

# Our right-hand-side function such that du/dt = f(u) for some vector
# u. Our vector u is (h,v,mp).
def rhs(u,t):
    """
    Returns the right-hand-side of the differential equation
    du/dt = f(u,t)
    for the rocket equations of motion
    """
    h = u[0]  # This is more verboes, but I think it's clearer.
    v = u[1]
    mp = u[2]
    hprime = v # dh/dt
    vprime = (-(ms+mp)*g + burn_rate(t)*ve - 0.5*rho*v*abs(v)*A*CD)/(ms+mp)
    mpprime = -burn_rate(t)
    return numpy.array([hprime,vprime,mpprime])


# We just steal the Euler Step method from the phugoid notebook.
# Only now it depends on time too!
def euler_step(u, t, f, dt):
    """
    Returns the solution at the next time-step using Euler's
    method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    t  : float
        the time at the previous time-step.
    f : function
        function to compute the right hand-side of the system
        of equation.
    dt : float
        time-increment.
    
    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """
    
    return u + dt * f(u,t)


# A method to perform an Euler solution.
def solve_rocket_equations():
    """
    Solves the rocket equations given an initial height h0
    and an initial velocity v0.

    Uses a first-order forward Euler method with a time step of dt
    and evolves from time t = 0 until the rocket crashes

    Returns a pair of arrays, the times t and the solutions u(t)
    """
    # initial data
    u0 = numpy.array([h0,v0,mp0])
    t = [0.0] # The array of times
    u = [u0]  # The array of evolved solutions
    # Evolve!
    while u[-1][0] >= 0.0:
        t.append(t[-1]+dt)
        u.append(euler_step(u[-1],t[-1],rhs,dt))
    return numpy.array(t),numpy.array(u)

def main():
    t,u=solve_rocket_equations()
    velocities = u[...,1]
    altitudes = u[...,0]
    max_velocity = numpy.max(velocities)
    mv_index = list(velocities).index(max_velocity)
    mv_time = t[mv_index]
    mv_h = u[mv_index,0]
    max_altitude = numpy.max(altitudes)
    ma_index = list(altitudes).index(max_altitude)
    ma_time = t[ma_index]
    final_time = t[-1]
    final_velocity = u[-1,1]
    print "We have {} time steps.".format(len(t))
    print "At t = 3.2s, the mass in the rocket is: {}".format(u[32,2])
    print "The maximum speed of the rocket is: {}".format(max_velocity)
    print "This occurs at time t = {} s".format(mv_time)
    print "The altitude at this time is h = {} m".format(mv_h)
    print "The rocket's maximum altitude is h = {} m".format(max_altitude)
    print "This occurs at time t = {} s.".format(ma_time)
    print "The rocket hits the ground at t = {} s.".format(final_time)
    print "The velocity when it hits the ground is v = {} m/s.".format(final_velocity)


if __name__ == "__main__":
    main()

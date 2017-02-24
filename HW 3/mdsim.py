#!/usr/bin/env python
from __future__ import print_function
import numpy as np

# ------------------------------------------------------------------------
# PHY466/MSE485 Atomic Scale Simulations
# Homework 2: Introduction to Molecular Dynamics
# ------------------------------------------------------------------------

from particleset import ParticleSet
"""
The ParticleSet class is designed to hold the position, velocity and accelerations of a set of particles. Initialization methods init_pos_cubic(cube_length) and init_vel(temperature) are provided for your convenience.

pset = ParticleSet(natom) will initialize a particle set pset
pset.size() will return the number of particles in pset

| --------------------------------------------------------- | ------------------- |
|                   for access of                           |    use method       |
| --------------------------------------------------------- | ------------------- |
| all particle positions in an array of shape (natom,ndim)  |    pset.all_pos()   |
| all particle velocities                                   |    pset.all_vel()   |
| all particle accelerations                                |    pset.all_accel() |
| particle i position in an array of shape (ndim)           |    pset.pos(i)      |
| particle i velocity                                       |    pset.vel(i)      |
| particle i acceleration                                   |    pset.accel(i)    |
| --------------------------------------------------------- | ------------------- |

| ----------------------------- | ------------------------------------ | 
|           to change           |             use method               |
| ----------------------------- | ------------------------------------ |
| all particle positions        |  pset.change_all_pos(new_pos_array)  |
| particle i position           |  pset.change_pos(i,new_pos)          |
| ditto for vel and accel       |  pset.change_*(i,new_*)              |
| ----------------------------- | ------------------------------------ |
"""

# Routines to ensure periodic boundary conditions that YOU must write.
# ------------------------------------------------------------------------
def pos_in_box(mypos, box_length):
    """ Return position mypos in simulation box using periodic boundary conditions. The simulation box is a cube of size box_length^ndim centered around the origin vec{0}. """
     
    new_pos = mypos.copy()
    ndim    = mypos.shape[0]
    # calculate new_pos
    for i in range(ndim):
        if new_pos[i] > box_length/2:
            new_pos[i]=new_pos[i]%(-box_length/2)
        elif new_pos[i] < (-box_length/2):
            new_pos[i]=new_pos[i]%(box_length/2)
    return new_pos
# def pos_in_box

def displacement(iat, jat, pset, box_length):
    """ Return the displacement of the iat th particle relative to the jat th particle. Unlike the distance function, here you will return a vector instead of a scalar """
    posi = pset.pos(iat)
    posj = pset.pos(jat)
    disp = posi.copy()
    ndim = posi.shape[0]
    for i in range(ndim):
        d = abs(posi[i]-posj[i])
        disp[i] = min(d,box_length-d)

    # calculate displacement of the iat th particle relative to the jat th particle
    # i.e. r_i - r_j
    # be careful about minimum image convention! i.e. periodic boundary
    return disp
# end def distance

def distance(iat, jat, pset, box_length):
    """ return the distance between particle i and j according to the minimum image convention. """

    dist = 0.0
    posi = pset.pos(iat)
    posj = pset.pos(jat)
    if posi.shape[0] != posj.shape[0]:
        print("Incorrect dimensions!")
        return
    ndim = posi.shape[0]
    A = np.zeros(ndim)
    for i in range(ndim):
        d = abs(posi[i]-posj[i])
        A[i] = min(d,box_length-d)

    dist = np.linalg.norm(A)
    # A = np.zeros(ndim)
    # B = np.zeros(ndim)
    # ndim = posi.shape[0]
    # #print(posi, posj)
    # # calculate distance with minimum image convention
    # # np.linalg.norm() may be useful here
    # for k in range(ndim):
    #     if posi[k]!=posj[k]:
    #         lmbda1 = (box_length/2 - posi[k])/(posi[k] - posj[k])
    #         lmbda2 = -(box_length/2 + posi[k])/(posi[k] - posj[k])
    #         for m in range(ndim):
    #             if k==m:
    #                 continue
    #             A[m] = lmbda1*(posi[m] - posj[m]) + posi[m]
    #             B[m] = lmbda2*(posi[m] - posj[m]) + posi[m]
    #             if A[m] > box_length/2:
    #                 A[m]=0
    #             if B[m] > box_length/2:
    #                 B[m] = 0

    return dist

# end def distance

# The Verlet time-stepping algorithm that YOU must write, dt is time step
# ------------------------------------------------------------------------
def verlet_next_pos(pos_t,vel_t,accel_t,dt):
    """
    We want to return position of the particle at the next moment t_plus_dt
    based on its position, velocity and acceleration at time t.  
    """
    pos_t_plus_dt = pos_t + vel_t*dt + 0.5*accel_t*dt**2

    return pos_t_plus_dt
# end def verlet_next_pos

def verlet_next_vel(vel_t,accel_t,accel_t_plus_dt,dt):
    """
    We want to return velocity of the particle at the next moment t_plus_dt, 
    based on its velocity at time t, and its acceleration at time t and t_plus_dt
    """
    vel_t_plus_dt = vel_t + 0.5*(accel_t+accel_t_plus_dt)*dt

    return vel_t_plus_dt
# end def verlet_next_vel

# We want Lennard-Jones forces. YOU must write this.
# ------------------------------------------------------------------------
def internal_force(iat,pset,box_length):
    """
    We want to return the force on atom 'iat' when we are given a list of 
    all atom positions. Note, pos is the position vector of the 
    1st atom and pos[0][0] is the x coordinate of the 1st atom. It may
    be convenient to use the 'displacement' function above. For example,
    disp = displacement( 0, 1, pset, box_length ) would give the position
    of the 1st atom relative to the 2nd, and disp[0] would then be the x coordinate
    of this displacement. Use the Lennard-Jones pair interaction. Be sure to avoid 
    computing the force of an atom on itself.
    """

    pos = pset.all_pos()  # positions of all particles
    mypos = pset.pos(iat) # position of the iat th particle
    npart =  pos.shape[0]
    #print npart
    force = np.zeros(pset.ndim())
    # calculate force
    for a in range(npart):
        r = distance(iat,a,pset,box_length)
        if r == 0:
            continue
        disp = displacement(iat,a,pset,box_length)
        for i in range(pset.ndim()):
            force[i] = 24*disp[i]*(2/(r**7) - 1/r**4)

    return force
# end def internal_force

def compute_energy(pset,box_length):
    natom = pset.size()  # number of particles
    pos = pset.all_pos() # all particle positions 
    vel = pset.all_vel() # all particle velocies

    tot_kinetic   = 0.0
    tot_potential = 0.0 
    tot_kinetic = 0.5*pset.mass()*(vel**2).sum()
    #mypos = pos[0]
    for a in range(1,len(pos)):
        r = distance(0,a,pset,box_length)
        if r == 0:
            continue
        tot_potential = tot_potential + 4*(1/r**12 - 1/r**6)

    tot_energy = tot_kinetic + tot_potential
    return (tot_kinetic, tot_potential, tot_energy)
# end def compute_energy

if __name__ == '__main__':

    num_atoms   = 64
    mass        = 48.0
    temperature = 0.728
    box_length  = 4.2323167
    nsteps      = 10
    dt          = 0.01

    # create and initialize a set of particles
    pset = ParticleSet(num_atoms,mass)
    pset.init_pos_cubic(box_length)
    pset.init_vel(temperature)

    # molecular dynamics simulation loop
    for istep in range(nsteps):

        # calculate properties of the particles
        print(istep, compute_energy(pset,box_length))

        # update positions
        for iat in range(num_atoms):
            my_next_pos = verlet_next_pos( pset.pos(iat), pset.vel(iat), pset.accel(iat), dt)
            new_pos = pos_in_box(my_next_pos,box_length)
            pset.change_pos(iat,new_pos)
        # end for iat
        
        # Q/ When should forces be updated?
        new_accel = pset.all_accel()

        # update velocities
        for iat in range(num_atoms):
            my_next_vel = verlet_next_vel( pset.vel(iat), pset.accel(iat), new_accel[iat], dt )
            pset.change_vel( iat, my_next_vel )
        # end for iat

    # end for istep

# end __main__

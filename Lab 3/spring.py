#!/usr/bin/env python

# Chris Dudiak
# This program investigates numerically the motion of a spring.

import math
import numpy as np
import operator as op
import sys
import matplotlib.pyplot as plt

# global constants
N = 4000

# explicit Euler method
def spring(show, x0, v0):
    h = .01
    # create a list of the two lists zipped together
    z = zip([x0 for i in range(N)], [v0 for i in range(N)])
    z0 = z[:]
    
    springProp(z, h)

    # unzip the lists back to x, v
    (xs, vs) = zip(*z)
    # plot the values
    fig = plt.figure()
    ts = np.arange(0, N) * h
    plt.plot(ts, xs, label="Position")
    plt.plot(ts, vs, label="Velocity")
    title = "Position and Velocity of Spring Over Time"
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Function Value")
    plt.legend()
    if show:
        plt.show()
    fig.savefig("explicit.png")
    
    # analytic
    xE = 5.0 * np.cos(ts)
    vE = -5.0 * np.sin(ts)
    # plot errors
    fig = plt.figure()
    plt.plot(ts, xE - xs, label="Position Error")
    plt.plot(ts, vE - vs, label="Velocity Error")
    title = "Position and Velocity Global Errors"
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Difference Between Analytic and Calculated")
    plt.legend()
    if show:
        plt.show()
    fig.savefig("errors.png")
    
    h2 = .0027
    
    # truncation error
    hs = [h2 * (0.5 ** i) for i in range(5)]
    errors = []
    for h in hs:
        xsForH = map(lambda x: x[0], springProp(z0[:], h))
        tvals = np.arange(0, N) * h
        xEtest = 5.0 * np.cos(tvals)
        errors.append(max(xEtest - xsForH))
    # plot errors
    fig = plt.figure()
    plt.plot(hs, errors)
    title = "Truncations Errors in Position Relative to h"
    plt.title(title)
    plt.xlabel("h")
    plt.ylabel("Maximum Error from Analytic Result")
    if show:
        plt.show()
    fig.savefig("trunc.png")
    
    # energy
    energy = map(lambda x, v: x**2 + v**2, xs, vs)
    # plot errors
    fig = plt.figure()
    plt.plot(ts, energy)
    title = "Normalized Total Energy vs Time"
    plt.title(title)
    plt.xlabel("Time, t")
    plt.ylabel("Normalized Total Energy")
    if show:
        plt.show()
    fig.savefig("energy.png")
    
    # implicit
    val = np.matrix([x0, v0]).transpose()
    h = .01
    mat = np.linalg.inv(np.matrix([[1, -h],[h, 1]]))
    impVals = [0 for i in range(N)]
    for i in range(N):
        impVals[i] = val
        val = mat * val
    xsImp = map(lambda x: x[0, 0], impVals)
    vsImp = map(lambda v: v[1, 0], impVals)
    
    # plot explicit phase
    fig = plt.figure()
    plt.plot(xs, vs)
    title = "Explicit Phase Space"
    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Velocity")
    if show:
        plt.show()
    fig.savefig("ephase.png")
    
    # plot implicit phase
    fig = plt.figure()
    plt.plot(xsImp, vsImp)
    title = "Implicit Phase Space"
    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Velocity")
    if show:
        plt.show()
    fig.savefig("iphase.png")
    
    # show difference with explicit
    # plot errors
    fig = plt.figure()
    plt.plot(ts, xE - xs, label="Position Error Explicit")
    plt.plot(ts, xE - xsImp, label="Position Error Implicit")
    title = "Position Global Errors"
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Difference Between Analytic and Calculated")
    plt.legend()
    if show:
        plt.show()
    fig.savefig("positionComp.png")
    
    fig = plt.figure()
    plt.plot(ts, vE - vs, label="Velocity Error Explicit")
    plt.plot(ts, vE - vsImp, label="Velocity Error Implicit")
    title = "Velocity Global Errors"
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Difference Between Analytic and Calculated")
    plt.legend()
    if show:
        plt.show()
    fig.savefig("velocityComp.png")
    
    energyImp = map(lambda x, v: x**2 + v**2, xsImp, vsImp)
    fig = plt.figure()
    plt.plot(ts, energy, label="Energy Explicit")
    plt.plot(ts, energyImp, label="Energy Implicit")
    title = "Velocity Global Errors"
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Difference Between Analytic and Calculated")
    plt.legend()
    if show:
        plt.show()
    fig.savefig("energyComp.png")
    
    # symplectic Euler
    h = .5
    ts = np.arange(0, N) * h
    z = z0[:]
    springSymp(z, h)
    # unzip the lists back to x, v
    (xsSym, vsSym) = zip(*z)
    # analytic
    xE = 5.0 * np.cos(ts)
    vE = -5.0 * np.sin(ts)
    # plot the values alongside analytic
    fig = plt.figure()
    plt.plot(xsSym, vsSym, label="Symplectic Phase")
    plt.plot(xE, vE, label="Exact Phase")
    title = "Symplectic Phase Space"
    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Velocity")
    plt.legend()
    if show:
        plt.show()
    fig.savefig("sphase.png")
    
    # symplectic Euler energy
    h = .01
    ts = np.arange(0, N) * h
    z = z0[:]
    springSymp(z, h)
    # unzip the lists back to x, v
    (xsSym, vsSym) = zip(*z)
    energySymp = map(lambda x, v: x**2 + v**2, xsSym, vsSym)
    # plot energy
    fig = plt.figure()
    plt.plot(ts, energySymp)
    title = "Normalized Total Energy vs Time"
    plt.title(title)
    plt.xlabel("Time, t")
    plt.ylabel("Normalized Total Energy")
    if show:
        plt.show()
    fig.savefig("symenergy.png")
    
    
# Takes a zipped list z and a step h and calculates all the new values for the system.
def springProp(z, h):
    # loop over the lists simultaneously and compute new values
    for i, (x, v) in enumerate(z):
        if (i + 1 >= N):
            break
        z[i + 1] = (x + h * v, v - h * x)
    return z
    
# Takes a zipped list z and a step h and calculates all the new values for the system
# using the symplectic Euler method.
def springSymp(z, h):
    # loop over the lists simultaneously and compute new values
    for i, (x, v) in enumerate(z):
        if (i + 1 >= N):
            break
        x_next = x + h * v
        z[i + 1] = (x_next, v - h * x_next)
    return z
        
    
if __name__ == "__main__":
    if len(sys.argv) != 4 and len(sys.argv) != 2:
        print "usage: %s bool_to_show_plot(0 or 1) [x0 v0]" % sys.argv[0]
        sys.exit(1)
    elif len(sys.argv) == 4:
        spring(int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]))
    else:
        spring(int(sys.argv[1]), 5.0, 0.0)
#!/usr/bin/env python

# Chris Dudiak
# This program contains various methods to solve for the roots of functions.

import math
import numpy as np
import operator as op
import sys
import matplotlib.pyplot as plt

# global vars
precision = .00001

def roots(show):
    (zb, sb) = bisection(func, -math.pi / 2.0, math.pi / 2.0)
    (zn, sn) = newtonRaphson(func, funcP, 0)
    (zs, ss) = secant(func, -math.pi / 2.0, math.pi / 2.0)
    print "Bisection Guess: %f" % zb
    print "Newton-Rpahson Guess: %f" % zn
    print "Secant Guess: %f" % zs
    
    xB = range(len(sb))
    xN = range(len(sn))
    xS = range(len(ss))
    
    yB = map(lambda y: math.log(abs(y)), sb)
    yN = map(lambda y: math.log(abs(y)), sn)
    yS = map(lambda y: math.log(abs(y)), ss)
    
    # plot the convergence on a log-log plot
    fig = plt.figure()
    plt.plot(xB, yB, label="Bisection")
    plt.plot(xN, yN, label="Newton-Raphson")
    plt.plot(xS, yS, label="Secant")
    title = "Convergence of Various Methods on sin(x) - .76"
    plt.title(title)
    plt.xlabel("Steps")
    plt.ylabel("Log f(x)")
    plt.legend()
    if show:
        plt.show()
    fig.savefig("convergence.png")
    
    # Question 3 - Orbit
    e = 0.617139
    T = 27906.98161
    c = 299792458
    a = 2.34186 * c
    
    orbP = lambda x : T / (2 * math.pi) * (1 - e * np.cos(x))
    xs = []
    ys = []
    deltaT = .001
    xsDelta = []
    ysDelta = []
    ts = range(-14000, 28000, 100)
    for t in ts:
        orb = lambda x : T / (2 * math.pi) * (x - e * np.sin(x)) - t
        orb2 = lambda x : T / (2 * math.pi) * (x - e * np.sin(x)) - (t + deltaT)
        (zero, _) = newtonRaphson(orb, orbP, t)
        (zero2, _) = newtonRaphson(orb2, orbP, t + deltaT)
        x = a * (np.cos(zero) - e)
        y = a * math.sqrt(1 - e**2) * np.sin(zero)
        xs.append(x)
        ys.append(y)
        
        xD = a * (np.cos(zero2) - e)
        yD = a * math.sqrt(1 - e**2) * np.sin(zero2)
        xsDelta.append(xD)
        ysDelta.append(yD)
        
    # plot the orbit
    fig = plt.figure()
    plt.plot(xs, ys, label="Orbit")
    title = "Orbit of Binary Pulsar 1913+16"
    plt.title(title)
    plt.xlabel("X position")
    plt.ylabel("Y Position")
    plt.legend()
    if show:
        plt.show()
    fig.savefig("orbit.png")
    
    # Part 4 - Velocity Curve
    th = -math.pi / 2.0
    rads = []
    for i in range(len(xs)):
        xp = (xsDelta[i] - xs[i]) / deltaT
        yp = (ysDelta[i] - ys[i]) / deltaT
        radV = np.dot([xp, yp], [np.cos(th), np.sin(th)]) / 1000.0
        rads.append(radV)
    
    tTs = map(lambda t : float(t) / T, ts)
    # plot the velocity curve
    fig = plt.figure()
    plt.plot(tTs, rads)
    title = "Velocity Curve of Binary Pulsar 1913+16"
    plt.title(title)
    plt.xlabel("Phase")
    plt.ylabel("Radial Velocity, km/s")
    if show:
        plt.show()
    fig.savefig("velocity.png")
       
    print "Found a Qualitative Velocity Match at theta = %f" % th    
    print "Successfully Ran Script"

def func(x):
    return np.sin(x) - .76
    
def funcP(x):
    return np.cos(x)

# Takes a function and bounds x1,x2 to locate the zero by bisection.
# returns the tuple (x0, [steps taken]) 
def bisection(f, x1, x2):
    acc = abs(x1 - x2)
    steps = []
    x0 = (x1 + x2) / 2.0
    while acc > precision:
        x0 = (x1 + x2) / 2.0
        guess = f(x0)
        if sign(guess) == sign(f(x1)):
            x1 = x0
        else:
            x2 = x0
        acc = abs(x1 - x2)
        steps.append(guess)
    return (x0, steps)

# Takes a function, its derivative, and initial guess x1 to locate 
# the zero by Newton-Raphson.
# returns the tuple (x0, [steps taken]) 
def newtonRaphson(f, fp, x1):
    guess = f(x1)
    steps = [guess]
    x2 = x1
    while abs(guess) > precision:
        x2 = x1 - guess / fp(x1)
        guess = f(x2)
        x1 = x2
        steps.append(guess)
    return (x2, steps)
    
# Takes a function, two initial guesses x1,x2 to locate 
# the zero by Secant method.
# returns the tuple (x0, [steps taken])
def secant(f, x1, x2):
    guess = f(x2)
    prev = f(x1)
    steps = [guess]
    while abs(guess) > precision:
        next = x2 - guess * (x2 - x1) / (guess - prev)
        prev = guess
        guess = f(next)
        steps.append(guess)
        x1 = x2
        x2 = next
    return (x2, steps)

# determines the sign of x
def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0           
    
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "usage: %s bool_to_show_plot(0 or 1) " % sys.argv[0]
        sys.exit(1)
    print "Running Methods on default function: f(x) = sin(x) - .76"
    roots(int(sys.argv[1]))
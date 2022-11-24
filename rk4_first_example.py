from rungeKutta import *

"""
Find the solution for the second order differential equation
y'' + 4y - cos(3x) = 0
x0 = 0
with initial conditions y(x0) = y(0) = 1 and y'(x0) = y'(0) = 2
using The Fourth Order-Runge Kutta Method
this works by splitting the problem into 2 first order differential equations:

y' = z = f(x, y, z)
z' = g(x, y, z)
with initial conditions y(x0) = y(0) = 1 and y'(0) = z(0) = 2

"""

if __name__ == '__main__':
    # set initial conditions y(x0) = y0, z(x0) = z0 and other parametrs
    # xn is the final point
    # h1 > 0 is step-size for Numerical Solution
    # h2 > 0 is step-size for Analytical Solution
    x0 = 0
    xn = 7
    h1 = 0.2
    h2 = 0.001
    y0 = 1
    z0 = 2

    xlist1, ylist1 = RungeKutta4(x0, xn, y0, z0, h1)

    xlist2 = [x*h2 for x in range(x0, int(xn/h2+1))]
    ylist2 = [exact(x) for x in xlist2]

    plot(xlist1, ylist1, xlist2, ylist2, 'blue', 'red', 'solid', 'solid', h1)

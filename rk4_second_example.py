"""
Find the solution for the second order differential equation
y'' - 4y = 0
x0 = 0
with initial conditions y(x0) = y(0) = 1 and y'(x0) = y'(0) = 2
using The Fourth Order-Runge Kutta Method
this works by splitting the problem into 2 first order differential equations:

y' = z = f(x, y, z)
z' = g(x, y, z)
with initial conditions y(x0) = y(0) = 1 and y'(0) = z(0) = 2

"""

from math import sqrt, exp, cos, sin
import matplotlib.pyplot as plt


def RungeKutta4(x0, xn, y0, z0, h):

    n = int((xn - x0)/h)
    # Containers for solutions
    xlist = [0] * (n + 1)
    ylist = [0] * (n + 1)
    zlist = [0] * (n + 1)

    xlist[0] = x = x0
    ylist[0] = y = y0
    zlist[0] = z = z0

    for i in range(1, n + 1):
        # see https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
        k1 = h * f(x, y, z)
        l1 = h * g(x, y, z)
        k2 = h * f(x + 0.5 * h, y + 0.5 * k1, z + 0.5*l1)
        l2 = h * g(x + 0.5 * h, y + 0.5 * k1, z + 0.5*l1)
        k3 = h * f(x + 0.5 * h, y + 0.5 * k2, z + 0.5*l2)
        l3 = h * g(x + 0.5 * h, y + 0.5 * k2, z + 0.5*l2)
        k4 = h * f(x + h, y + k2, z + l2)
        l4 = h * g(x + h, y + k2, z + l2)
        xlist[i] = x = x0 + i * h
        ylist[i] = y = y + (k1 + 2*k2 + 2*k3 + k4) / 6
        zlist[i] = z = z + (l1 + 2*l2 + 2*l3 + l4) / 6

    return xlist, ylist


def g(x, y, z):
    # z' = 4y
    return 4*y


def f(x, y, z):
    # y' = z
    return z


def exact(x):
    # analytical solution y = exp(2*x)
    return exp(2*x)


def plot(x1, y1, x2, y2, color1, color2, linestyle1, linestyle2, h):

    dpi = 80
    fig = plt.figure(dpi=dpi, figsize=(1600 / dpi, 900 / dpi))

    plt.plot(x1, y1, color=color1, linestyle=linestyle1,
             label='Numerical Solution')
    plt.plot(x2, y2, color=color2, linestyle=linestyle2,
             label='Analytical Solution')
    equation = "y'' - 4y = 0"
    name = "Gordievich Kirill"
    university = "SMTU"
    title = '%s, h = %s, The Runge–Kutta method by %s, %s' %(equation, h, name, university)
    plt.title(title, fontsize=20)
    plt.xlabel('x', fontsize=20)
    plt.ylabel('y', fontsize=20)
    plt.legend(loc='upper right')
    plt.show()
    # Uncomment the following to print the figure:
    #file_name = str(h1) + '_2_RungeKutta4.png'
    #fig.savefig(file_name)

if __name__ == '__main__':
    # set initial conditions y(x0) = y0, z(x0) = z0 and other parametrs
    # xn is the final point
    # h1 > 0 is step-size for Numerical Solution
    # h2 > 0 is step-size for Analytical Solution
    x0 = 0
    xn = 3
    h1 = 0.5
    h2 = 0.001
    y0 = 1
    z0 = 2

    xlist1, ylist1 = RungeKutta4(x0, xn, y0, z0, h1)

    xlist2 = [x*h2 for x in range(x0, int(xn/h2+1))]
    ylist2 = [exact(x) for x in xlist2]

    plot(xlist1, ylist1, xlist2, ylist2, 'blue', 'red', 'solid', 'solid', h1)


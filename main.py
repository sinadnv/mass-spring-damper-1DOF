import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

# Define the variables
m = 1
c = 1
k = 4
F0 = 10
x0 = [0, 0]


#Define the function of interest.
# I would like to plot the response of the system as a function of the frequency. so I need to include omega as
# a function input as well.
def mass_spring_damper(x,t,omega):
    F = F0 * np.sin(omega* t)
    x1dot = x[1]
    x2dot = 1/m*(F-c*x[1]-k*x[0])
    return [x1dot, x2dot]

#Define the range of analysis
tmax = 50
omegamax=10
dt = .005
# dt is an extremely important factor to determine the magnitude the vibrations.
# It's best to check the energy of the system in the discretization method to check for convergence
t = np.linspace(0,tmax,int(tmax/dt))
omega = np.linspace(.1,omegamax,omegamax*10)

# Solve ODE for all excitation frequencies, w
umax = []
for w in omega:
    sol = odeint(mass_spring_damper, x0, t, args=(w,))
    u = sol[:, 0]
    umax.append(max(u))

#plt.plot(omega,umax)

x1 = np.zeros(len(t))
x2 = np.zeros(len(t))
u_discr = []
for w in omega:
    F = F0*np.sin(w*t)
    for i in range(len(t)-1):
        x1[i+1] = x2[i]*dt+x1[i]
        x2[i+1] = dt/m*F[i]+(1-c*dt/m)*x2[i]-k*dt/m*x1[i]
    u_discr.append(max(x1))

plt.plot(omega, u_discr)
plt.xlabel('$\omega$')
plt.ylabel('A')
plt.title('Frequency Response of 1-DOF system')
plt.grid()
plt.show()

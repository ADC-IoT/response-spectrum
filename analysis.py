import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# --initialization
# M : ton
# C : KNsec/m
# K : KN/m
# dt: sec
# wave : cm/sec2

# Parameter setting
# Analysis
beta = 1/4

# Frame
h = 0.03
m1 = 100

# System
N = 1

# Wave
Ndiv = 1
dt = 0.02
ddt = dt/float(Ndiv)
Ndata = 2688
EQ = 0.01  # cm to m

dis1 = 0
vel1 = 0
acc1 = 0

# INPUT FILE
file = open('ELC-NS.csv', 'r', encoding="utf-8")
A0 = np.loadtxt(file, usecols=(0,), delimiter=',', skiprows=3)
F = -m1*A0

# list for data plotting
xaxis = []
wave = []
rac=np.zeros(Ndata)
vel=np.zeros(Ndata)
dis=np.zeros(Ndata)
Amax = []
Vmax = []
Xmax = []
Tr = []

def nmk(beta,m1,c1,k1,dt,ff,acc1,vel1,dis1):
    acc2=(ff-c1*(vel1+0.5*dt*acc1)-k1*(dis1+dt*vel1+(0.5-beta)*dt*dt*acc1))/(m1+c1*0.5*dt+k1*beta*dt*dt)
    vel2=vel1+0.5*dt*(acc1+acc2)
    dis2=dis1+dt*vel1+(0.5-beta)*dt*dt*acc1+beta*dt*dt*acc2
    return acc2,vel2,dis2


# calculation
for t in range(1,500):
    T = 0.01 * t
    Tr.append(T)
    k1 = 4*np.pi**2*m1/T**2
    c1 = 2*h*np.sqrt(k1*m1)

    for i in range(0,Ndata):
        f1=0.0
        if 1<=i: f1=F[i-1]
        f2=F[i]
        for k in range(0,Ndiv):
            fm=f1+(f2-f1)/float(Ndiv)*float(k)
            acc2,vel2,dis2=nmk(beta,m1,c1,k1,ddt,fm,acc1,vel1,dis1)
            acc1=acc2
            vel1=vel2
            dis1=dis2
        rac[i]=acc1+A0[i]
        vel[i]=vel1
        dis[i]=dis1

    Amax.append(np.max(np.abs(rac)))
    Vmax.append(np.max(np.abs(vel)))
    Xmax.append(np.max(np.abs(dis)))

# data plotting
plt.subplot(3,1,1)
plt.title("spectrum")
plt.xlabel('period (s)')
plt.ylabel('acc. (cm/sec2)')
plt.plot(Tr, Amax)

plt.subplot(3,1,2)
plt.xlabel('period (s)')
plt.ylabel('vel. (cm/sec)')
plt.plot(Tr, Vmax)

plt.subplot(3,1,3)
plt.xlabel('period (s)')
plt.ylabel('dis. (cm)')
plt.plot(Tr, Xmax)

plt.show()

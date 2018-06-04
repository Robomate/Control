#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Purpose: LQR Controller, RLC Series Circuit Simulation
Description: see 2017_10_24_01.m
Version: 5/2018 Roboball (MattK.)

Example Sketch: RLC Series Circuit

          L		 C		
 +   o---###----###-----------o 
     |                 |      |
     |                 #      |
    V_in               # R  V_out
     |                 #      |
     |                 |      |
GND  o                 o      o 

"""
import matplotlib.pyplot as plt
import math as m
import numpy as np
from numpy import linalg as la
from control.matlab import *    # MATLAB-like functions
import control as pc


# https://github.com/python-control/python-control/tree/master/examples

########################
# Init Globals
########################

# System parameters
L = 1e-3 # Vs/A
C = 1e-6 # As/V
R = 10   # V/A

##############################
# Tf Function Representation
##############################
sys2tf = tf([1], [L * C, R * C,1])
print(sys2tf)

# Plot step response
(y1a, T1a) = step(sys2tf)              # not correct!!
fig1 = plt.figure(1, figsize=(8,6))
plt.title('Step Response')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (V)")
plt.plot(T1a, y1a, color='#20B2AA',label='TF')  # why 1e9
#plt.axis([0,len(T1a),0,max(y1a)]) # opt.: set ticks
plt.xlim(xmin=0) 
plt.ylim(ymin=0) 
plt.grid(True)
plt.hold(True)

# Plot bode plot
fig2 = plt.figure(2, figsize=(8,6))
mag, phase, omega = pc.bode_plot(sys2tf,dB= True,  Hz=False) # not correct!!

# tf to state space
ss_check = tf2ss(sys2tf)
print(ss_check)

##############################
# State Space Representation
##############################
A = [[-R/L,-1/L],
     [1/C ,0   ]]
b = [[1/L],
     [0  ]]
c = [[0,1]]
d = [[0  ]]

sys1ss = ss(A, b, c, d)
print(sys1ss)
fig2 = plt.figure(3, figsize=(8,6))
mag, phase, omega = pc.bode_plot(sys1ss,dB= True,  Hz=False) 
#print(sys1ss)
# Plot step response state space
(y2a, T2a) = step(sys1ss)
plt.figure(1, figsize=(8,6))
plt.plot(T2a, y2a, 'bs',label='SS')  # why 1e9

# Regelungsnormalform
An = [[0     ,1   ],
      [-1/L/C,-R/L]]
bn = [[0   ],
     [1/L/C]]      
cn = [[1,0]]
dn = [[0  ]]
sys2ss = ss(An, bn, cn, dn)
print(sys2ss)
(y3a, T3a) = step(sys2ss)
plt.figure(1, figsize=(8,6))
plt.plot(T3a, y3a, 'm--',label='SS NF-Form')

# Jordan Normalform
lambdas, M = la.eig(np.array(A)) # Eigen_vectors, Eigen_values

T = la.inv(M)
A_ = np.dot(T,np.dot(np.array(A),M))
b_ = np.dot(T,np.array(b))
c_ = np.dot(np.array(c),M)

tend = 1.4e-3;
T4a = np.arange(0,tend,tend/100)
u0 = 1;
x1 = b_[0]*u0/-lambdas[0]*(1-np.exp(lambdas[0]*T4a))
x2 = b_[1]*u0/-lambdas[1]*(1-np.exp(lambdas[1]*T4a))

y4a = np.real(c_[0][0] * x1 + c_[0][1] * x2)
#~ print(np.real(y4a))
plt.figure(1, figsize=(8,6))
plt.plot(T3a, y4a, 'r+',label='Jordan NF-Form')

# Discretisation
f0 = (np.real(lambdas[0])**2+np.imag(lambdas[1])**2)**0.5/2/m.pi
T0 = 1/f0
Ts = T0/10 # time step
#print(Ts)
Phi = np.eye(len(A)) 
gamma = np.eye(len(A))
for i in range(1,21):
	Phi = Phi + np.linalg.matrix_power(A,i)*Ts**i/m.factorial(i)
	gamma = gamma + np.linalg.matrix_power(A,i)*Ts**i/m.factorial(i+1)
gamma = Ts * np.dot(gamma,b)

# double check with c2d
sys1dss = c2d(sys1ss ,Ts)
print()
print(sys1dss)
T5a = np.arange(0,tend,Ts)
x = np. zeros_like(np.array(b))
y5a = np. zeros_like(T5a)

for j in range(1,len(T5a)+1):
	y5a[j-1] = np.dot(c,x)
	x = np.dot(Phi,x) + gamma * u0

plt.figure(1, figsize=(8,6))
plt.step(T5a, y5a, 'g',label='Discrete SS')


# put at the end to show plots
plt.legend(loc='upper right')
#plt.show() 







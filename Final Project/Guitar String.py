# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 20:18:58 2021

@author: brand
"""

import numpy as np
import matplotlib.pyplot as plt
from control.matlab import *

s = tf('s')
#E(plain), B(plain), G(plain), D(Wound), A(wound), E(wound)
r = np.array([0.2540, 0.3302, 0.4310, 0.6604, 0.9144, 1.1684])*10**-3 

E = 200E9 # Steel Young Modulus
rr=r[0]
rho = 8000 #Steel Density
R = 10

#Number of nodes:
#n = input("Type in an odd numbers of interior nodes: ");
n = 5;

#Type in the length of string:
#l = input("\nType in the length of the string: ");
l = 10;

#Derived Parameter
dl = l/(n+1);
A = np.pi*rr**2;
k = E*A/dl;
R = R/n;
m = rho*A*l/n;

#Location of the Interior and Boundary Nodes
loc = np.zeros((n+2,1));
loc[0] = 0;
loc[n+1] = l;

for ii in range(1,n+1):
    loc[ii] = ii*dl
    
#Z,Y,H T.F
#Interiors Nodes are evenly spaced
#Distance b/w boundary and interior nodes are 1/2 of the interior nodes'
Z = []
for ii in range(0,n+1):
    Z.append([])
    for jj in range(0,n+1):
        Z[ii].append(0*s)
        
        
if n == 1:
    Z = m*s+k/s+R;
else :
    #For 1st boundary node
    Z[0][0] = m*s + 2*k/s + 2*R;
    Z[0][1] = -(k/s + R);
    for ii in range(1,n-2):
        Z[ii][ii-1] = Z[ii-1][ii];
        Z[ii][ii] = m*s -2*Z[ii][ii-1];
        Z[ii][ii+1] = Z[ii][ii-1];
    #For last boundary node
    Z[n-1][n-2] = Z[1][2];
    Z[n-1][n-1] = Z[1][1]; 
Z = np.matrix(Z)
Y = np.linalg.inv(Z)
#H = Y/s;
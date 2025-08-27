from sage.all import *

from helper_functions import Rz

C1=vector([3939,0 ,4340])/5861
C2=vector([7855,4178,4484])/10**4
C3=vector([9526,2057,1102])/10**4

RUP= matrix(SR,90, 3) # fill a 90 x 3 matrix with the vertices in the order as described in the paper
for i in range(3):
    for k in range(15):
        for l in range(2):
            if i==0:
                RUP[k+15*i+45*l,:]=(-1)**l*Rz(2*pi*k/15)*C1
            if i==1:
                RUP[k+15*i+45*l,:]=(-1)**l*Rz(2*pi*k/15)*C2
            if i==2:
                RUP[k+15*i+45*l,:]=(-1)**l*Rz(2*pi*k/15)*C3

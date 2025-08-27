from sage.all import *

from helper_functions import Rz

C1=vector([152024884,0 ,210152163])/259375205
C2=vector([6632738028,6106948881,3980949609])/10**10
C3=vector([8193990033,5298215096,1230614493])/10**10

NOP = matrix(SR,90, 3) # fill a 90 x 3 matrix with the vertices in the order as described in the paper
for i in range(3):
    for k in range(15):
        for l in range(2):
            if i==0:
                NOP[k+15*i+45*l,:]=(-1)**l*Rz(2*pi*k/15)*C1
            if i==1:
                NOP[k+15*i+45*l,:]=(-1)**l*Rz(2*pi*k/15)*C2
            if i==2:
                NOP[k+15*i+45*l,:]=(-1)**l*Rz(2*pi*k/15)*C3



from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# This is for the fcc!

def dot(S1, S2):
    return S1[0]*S2[0]+S1[1]*S2[1]+S1[2]*S2[2]
    
def cross(S1, S2):
    #print "checking cross, [a,b,c], a=",S1[1]*S2[2]-S1[2]*S2[1], ", b=", S1[2]*S2[0]-S1[0]*S2[2], ", c =", S1[0]*S2[1]-S1[1]*S2[0]
    return [ S1[1]*S2[2]-S1[2]*S2[1], S1[2]*S2[0]-S1[0]*S2[2], S1[0]*S2[1]-S1[1]*S2[0]]
    
def normalize_vec(S):
    return 1./sqrt(S[0]*S[0]+S[1]*S[1]+S[2]*S[2])
    
def vecdiff(S1, S2):
    return [S2[0]-S1[0], S2[1]-S1[1], S2[2]-S1[2]]

def planecheck(P, Q, R, A, i):
    PQ = vecdiff(P,Q)
    QR = vecdiff(Q,R)
    n = cross(PQ, QR)
    eq = n[0]*(A[0]-P[0]) + n[1]*(A[1]-P[1]) + n[2]*(A[2]-P[2])
    if abs(eq)<1e-6:
        #print "Point ", i, " is in the plane!"
        return 1
    else:
        #print "Spin", i, "not in plane, eq = ", eq
        return 0


infile = open("fcc6x6x6_nnJxz1_nnJyz1_beta100_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_spinconfigsameprint.txt", "r")
periodic = 1 # 1: false, we have open BCS


firstline = infile.readline() # Reads the first line only # Should maybe expand this. But who knows
beta, dim, L1, L2, L3 = firstline.split()

beta = float(beta); dim = float(dim); L1 = int(L1); L2 = int(L2); L3 = int(L3);
if dim==1:
    N = L1
elif dim==2:
    N = L1*L2
elif dim==3:
    N = L1*L2*L3

print "beta = ", beta

# Getting lists ready to store the data
Sxes    = []
Sys     = []
Szs     = []

# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        spinput = float(words[3])
        Sxes.append(spinput) 
        spinput = float(words[4])
        Sys.append(spinput)
        spinput = float(words[5])
        Szs.append(spinput)
        
# Remember to close the file
infile.close()

# Converting to arrays
Sxes    = array(Sxes)
Sys     = array(Sys)
Szs     = array(Szs)

sxav = 0
syav = 0
szav = 0
xlargerthanz = 0
xequaltoz = 0
print "Checking the normalization:"
for i in range(len(Sxes)):
    sxav+=abs(Sxes[i])
    syav+=abs(Sys[i])
    szav+=abs(Szs[i])
    if(Sxes[i]>Szs[i]):
        xlargerthanz +=1
    if(Sxes[i]==Szs[i]):
        xequaltoz    +=1
    S = [Sxes[i], Sys[i], Szs[i]]
    #print normalize_vec(S)

sxav = sxav/(N)
szav = szav/(N)
syav = syav/(N)

print "No of times Sx was larger than Sz: ", xlargerthanz, " out of ", L1*L2*L3
print "No of times Sx was equal t0 Sz: ", xequaltoz, " out of ", L1*L2*L3

print "Average of abs(Sx): ", sxav
print "Average of abs(Sy): ", syav
print "Average of abs(Sz): ", szav


# Getting the neighbours in the xy-direction
infile = open("/home/ubu/Master/Spincorrelation_Results/Yneighbours/fcc6x6x6_xyneighbours_eachpoint.txt", "r")
# Getting lists ready to store the data
xyns  = []

# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        n = int(words[1])
        xyns.append(n)
        
xyns = array(xyns)

infile.close()

# Getting the neighbours in the xz-direction
infile = open("/home/ubu/Master/Spincorrelation_Results/Yneighbours/fcc6x6x6_xzneighbours_eachpoint.txt", "r")
# Getting lists ready to store the data
xzns  = []
point = []

# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        n = int(words[0])
        point.append(n)
        n = int(words[1])
        xzns.append(n)
        
xzns = array(xzns)
point = array(point)

infile.close()
    
# Getting the neighbours in the xy-direction
infile = open("/home/ubu/Master/Spincorrelation_Results/Yneighbours/fcc6x6x6_yzneighbours_eachpoint.txt", "r")
# Getting lists ready to store the data
yzns  = []

# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        n = int(words[1])
        yzns.append(n)
        
yzns = array(yzns)

infile.close()


# xy-neighbours
angleav = 0
angleav_abs = 0
dotav = 0
nxav = 0
nyav = 0
nzav = 0
mindot = 1000
maxdot = -1000
maxangle = -1000
minangle = 1000
dplessthan0 = []
ndplessthan0 = 0

myrange = len(xyns)
for i in range(myrange):
    n = point[i]
    nb = xyns[i]
    S1 = [Sxes[n], Sys[n], Szs[n]]
    S2 = [Sxes[nb], Sys[nb], Szs[nb]]
    dotres = dot(S1, S2)
    angle = arccos(dotres)
    angleav += angle
    dotav += dotres
    normvec = cross(S1, S2)
    nfactor = normalize_vec(normvec) # normalization factor
    # We want the normal vector to be of unit length
    normvec[0] = nfactor*normvec[0]
    normvec[1] = nfactor*normvec[1]
    normvec[2] = nfactor*normvec[2]      
    nxav += normvec[0]
    nyav += normvec[1]
    nzav += normvec[2]
    if dotres<mindot:
        mindot = dotres
    if dotres>maxdot:
        maxdot = dotres
    if angle<minangle:
        minangle = angle
    if angle>maxangle:
        maxangle = angle
    # Check if the dot product is less than zero
    if dotres<0:
        dplessthan0.append(dotres)
        ndplessthan0 += 1
    #print "Spin ", i, " and ", nb, ":"
    #print "Normvec", normvec
    #print "nfactor:", nfactor
    #print "The dot product between spin ", i, " and spin ", nb, ": ", dotres, ". Angle between them: ", angle
    #print "The cross product between spin ", i, " and spin ", nb, ": ", normvec
    #print "Normalization factor for cross product:", nfactor

angleav     = angleav/myrange
angleav_abs = angleav_abs/myrange
dotav       = dotav/myrange
nxav        = nxav/myrange
nyav        = nyav/myrange
nzav        = nzav/myrange
nlengthav   = sqrt(nxav*nxav+nyav*nyav+nzav*nzav)
print " "
print "Neighbours in the xy-direction:"
print "Average angle between adjacent spins: ", angleav
print "Delta theta/theta_av: ", (maxangle-minangle)/angleav
print "Average dot product between adjacent spins: ", dotav
print "Average normal vector to spin plane: [", nxav, ",", nyav, ",", nzav, "]"
print "Length of average normal vector:", nlengthav 
print "Minimal value of dot product:", mindot
print "Maximal value of dot product:", maxdot

print "beta = ", beta

print "The dot product is less than zero for:"
for i in range(ndplessthan0):
    print "dotprod i =", dplessthan0[i]
if ndplessthan0==0:
    print "None"


# xz-neighbours
angleav = 0
angleav_abs = 0
dotav = 0
nxav = 0
nyav = 0
nzav = 0
mindot = 1000
maxdot = -1000
maxangle = -1000
minangle = 1000
dplessthan0 = []
ndplessthan0 = 0

myrange = len(xzns)
for i in range(myrange):
    n = point[i]
    nb = xzns[i]
    S1 = [Sxes[n], Sys[n], Szs[n]]
    S2 = [Sxes[nb], Sys[nb], Szs[nb]]
    dotres = dot(S1, S2)
    angle = arccos(dotres)
    angleav += angle
    dotav += dotres
    normvec = cross(S1, S2)
    nfactor = normalize_vec(normvec) # normalization factor
    # We want the normal vector to be of unit length
    normvec[0] = nfactor*normvec[0]
    normvec[1] = nfactor*normvec[1]
    normvec[2] = nfactor*normvec[2]      
    nxav += normvec[0]
    nyav += normvec[1]
    nzav += normvec[2]
    if dotres<mindot:
        mindot = dotres
    if dotres>maxdot:
        maxdot = dotres
    if angle<minangle:
        minangle = angle
    if angle>maxangle:
        maxangle = angle
    # Check if the dot product is less than zero
    if dotres<0:
        dplessthan0.append(dotres)
        ndplessthan0 += 1
    #print "Spin ", i, " and ", nb, ":"
    #print "Normvec", normvec
    #print "nfactor:", nfactor
    #print "The dot product between spin ", i, " and spin ", nb, ": ", dotres, ". Angle between them: ", angle
    #print "The cross product between spin ", i, " and spin ", nb, ": ", normvec
    #print "Normalization factor for cross product:", nfactor

angleav     = angleav/myrange
angleav_abs = angleav_abs/myrange
dotav       = dotav/myrange
nxav        = nxav/myrange
nyav        = nyav/myrange
nzav        = nzav/myrange
nlengthav   = sqrt(nxav*nxav+nyav*nyav+nzav*nzav)
print "Neighbours in the xz-direction:"
print "Average angle between adjacent spins: ", angleav
print "Delta theta/theta_av: ", (maxangle-minangle)/angleav
print "Average dot product between adjacent spins: ", dotav
print "Average normal vector to spin plane: [", nxav, ",", nyav, ",", nzav, "]"
print "Length of average normal vector:", nlengthav 
print "Minimal value of dot product:", mindot
print "Maximal value of dot product:", maxdot

print "beta = ", beta

print "The dot product is less than zero for:"
for i in range(ndplessthan0):
    print "dotprod i =", dplessthan0[i]
if ndplessthan0==0:
    print "None"
    
# yz-neighbours
angleav = 0
angleav_abs = 0
dotav = 0
nxav = 0
nyav = 0
nzav = 0
mindot = 1000
maxdot = -1000
maxangle = -1000
minangle = 1000
dplessthan0 = []
ndplessthan0 = 0

myrange = len(yzns)
for i in range(myrange):
    n = point[i]
    nb = yzns[i]
    S1 = [Sxes[n], Sys[n], Szs[n]]
    S2 = [Sxes[nb], Sys[nb], Szs[nb]]
    dotres = dot(S1, S2)
    angle = arccos(dotres)
    angleav += angle
    dotav += dotres
    normvec = cross(S1, S2)
    nfactor = normalize_vec(normvec) # normalization factor
    # We want the normal vector to be of unit length
    normvec[0] = nfactor*normvec[0]
    normvec[1] = nfactor*normvec[1]
    normvec[2] = nfactor*normvec[2]      
    nxav += normvec[0]
    nyav += normvec[1]
    nzav += normvec[2]
    if dotres<mindot:
        mindot = dotres
    if dotres>maxdot:
        maxdot = dotres
    if angle<minangle:
        minangle = angle
    if angle>maxangle:
        maxangle = angle
    # Check if the dot product is less than zero
    if dotres<0:
        dplessthan0.append(dotres)
        ndplessthan0 += 1
    #print "Spin ", i, " and ", nb, ":"
    #print "Normvec", normvec
    #print "nfactor:", nfactor
    #print "The dot product between spin ", i, " and spin ", nb, ": ", dotres, ". Angle between them: ", angle
    #print "The cross product between spin ", i, " and spin ", nb, ": ", normvec
    #print "Normalization factor for cross product:", nfactor

angleav     = angleav/myrange
angleav_abs = angleav_abs/myrange
dotav       = dotav/myrange
nxav        = nxav/myrange
nyav        = nyav/myrange
nzav        = nzav/myrange
nlengthav   = sqrt(nxav*nxav+nyav*nyav+nzav*nzav)
print "Neighbours in the yz-direction:"
print "Average angle between adjacent spins: ", angleav
print "Delta theta/theta_av: ", (maxangle-minangle)/angleav
print "Average dot product between adjacent spins: ", dotav
print "Average normal vector to spin plane: [", nxav, ",", nyav, ",", nzav, "]"
print "Length of average normal vector:", nlengthav 
print "Minimal value of dot product:", mindot
print "Maximal value of dot product:", maxdot

print "beta = ", beta

print "The dot product is less than zero for:"
for i in range(ndplessthan0):
    print "dotprod i =", dplessthan0[i]
if ndplessthan0==0:
    print "None"

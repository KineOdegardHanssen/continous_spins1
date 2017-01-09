from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys


betas     = []          # List of beta values
# Beta = 100
# Open the file for reading
'''
infile = open("fcc10t10t10_iso1_beta100_acceptancerate.txt", "r")


# Getting lists ready to store the data
betas.append(100)
acceptancerates100 = []    # Acceptancerate for beta=100 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates100.append(ar)    
 
# Remember to close the file
infile.close()
'''

# Beta = 5
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta5_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(5)
acceptancerates5 = []    # Acceptancerate for beta=100 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates5.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 2.5
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta2p5_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(2.5)
acceptancerates2p5 = []    # Acceptancerate for beta=2.5 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates2p5.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 1
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta1_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(1)
acceptancerates1 = []    # Acceptancerate for beta=1 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates1.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 0.75
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta0p75_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(0.75)
acceptancerates0p75 = []    # Acceptancerate for beta=0.75 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates0p75.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 0.5
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta0p5_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(0.5)
acceptancerates0p5 = []    # Acceptancerate for beta=0.5 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates0p5.append(ar)    
 
# Remember to close the file
infile.close()


# Beta = 0.4
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta0p4_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(0.4)
acceptancerates0p4 = []    # Acceptancerate for beta=0.4 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates0p4.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 0.3
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta0p3_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(0.3)
acceptancerates0p3 = []    # Acceptancerate for beta=0.3 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates0p3.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 0.2
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta0p2_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(0.2)
acceptancerates0p2 = []    # Acceptancerate for beta=0.2 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates0p2.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 0.1
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta0p1_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(0.1)
acceptancerates0p1 = []    # Acceptancerate for beta=0.1 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates0p1.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 0.05
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta0p05_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(0.05)
acceptancerates0p05 = []    # Acceptancerate for beta=0.05 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates0p05.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 0.03
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta0p03_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(0.03)
acceptancerates0p03 = []    # Acceptancerate for beta=0.03 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates0p03.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 0.01
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta0p01_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(0.01)
acceptancerates0p01 = []    # Acceptancerate for beta=0.01 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates0p01.append(ar)    
 
# Remember to close the file
infile.close()

# Beta = 0
# Open the file for reading
infile = open("fcc10t10t10_iso1_beta0_acceptancerate.txt", "r")

# Getting lists ready to store the data
betas.append(0)
acceptancerates0 = []    # Acceptancerate for beta=0 #Is this too much? Well, could always comment it out
#                         # afterwards, anyway

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        ar = float(words[0])
        acceptancerates0.append(ar)    
 
# Remember to close the file
infile.close()

acceptancerates0     = array(acceptancerates0)
acceptancerates0p01  = array(acceptancerates0p01)
acceptancerates0p03  = array(acceptancerates0p03)
acceptancerates0p05  = array(acceptancerates0p05)
acceptancerates0p1   = array(acceptancerates0p1)
acceptancerates0p2   = array(acceptancerates0p2)
acceptancerates0p3   = array(acceptancerates0p3)
acceptancerates0p4   = array(acceptancerates0p4)
acceptancerates0p5   = array(acceptancerates0p5)
acceptancerates0p75  = array(acceptancerates0p75)
acceptancerates1     = array(acceptancerates1)
acceptancerates2p5   = array(acceptancerates2p5)
acceptancerates5     = array(acceptancerates5)
#acceptancerates100   = array(acceptancerates100)

N = len(acceptancerates1)

av0    = average(acceptancerates0)
av0p01 = average(acceptancerates0p01)
av0p03 = average(acceptancerates0p03)
av0p05 = average(acceptancerates0p05)
av0p1  = average(acceptancerates0p1)
av0p2  = average(acceptancerates0p2)
av0p3  = average(acceptancerates0p3)
av0p4  = average(acceptancerates0p4)
av0p5  = average(acceptancerates0p5)
av0p75 = average(acceptancerates0p75)
av1    = average(acceptancerates1) 
av2p5  = average(acceptancerates2p5) 
av5    = average(acceptancerates5) 
#av100  = average(acceptancerates1) 

std0    = 0
std0p01 = 0
std0p03 = 0
std0p05 = 0
std0p1  = 0
std0p2  = 0
std0p3  = 0
std0p4  = 0
std0p5  = 0
std0p75 = 0
std1    = 0 
std2p5  = 0
std5    = 0
#std100  = 0
for i in range(0, N):
    std0    += (acceptancerates0[i]-av0)*(acceptancerates0[i]-av0)
    std0p01 += (acceptancerates0p01[i]-av0p01)*(acceptancerates0p01[i]-av0p01)
    std0p03 += (acceptancerates0p03[i]-av0p03)*(acceptancerates0p03[i]-av0p03)
    std0p05 += (acceptancerates0p05[i]-av0p05)*(acceptancerates0p05[i]-av0p05)
    std0p1  += (acceptancerates0p1[i]-av0p1)*(acceptancerates0p1[i]-av0p1)
    std0p2  += (acceptancerates0p2[i]-av0p2)*(acceptancerates0p2[i]-av0p2)
    std0p3  += (acceptancerates0p3[i]-av0p3)*(acceptancerates0p3[i]-av0p3)
    std0p4  += (acceptancerates0p4[i]-av0p4)*(acceptancerates0p4[i]-av0p4)
    std0p5  += (acceptancerates0p5[i]-av0p5)*(acceptancerates0p5[i]-av0p5)
    std0p75 += (acceptancerates0p75[i]-av0p75)*(acceptancerates0p75[i]-av0p75)
    std1    += (acceptancerates1[i]-av1)*(acceptancerates1[i]-av1)
    std2p5  += (acceptancerates2p5[i]-av2p5)*(acceptancerates2p5[i]-av2p5)
    std5    += (acceptancerates5[i]-av5)*(acceptancerates5[i]-av5)
    #std100  += (acceptancerates100[i]-av100)*(acceptancerates100[i]-av100)
    
std0    = sqrt(1.0*std0/N)
std0p01 = sqrt(1.0*std0p01/N)
std0p03 = sqrt(1.0*std0p03/N)
std0p05 = sqrt(1.0*std0p05/N)
std0p1  = sqrt(1.0*std0p1/N)
std0p2  = sqrt(1.0*std0p2/N)
std0p3  = sqrt(1.0*std0p3/N)
std0p4  = sqrt(1.0*std0p4/N)
std0p5  = sqrt(1.0*std0p5/N)
std0p75 = sqrt(1.0*std0p75/N)
std1    = sqrt(1.0*std1/N)
std2p5  = sqrt(1.0*std2p5/N)
std5    = sqrt(1.0*std5/N)


acceptanceratesvsbeta = [av5, av2p5, av1, av0p75, av0p5, av0p4, av0p3, av0p2, av0p1, av0p05, av0p03, av0p01, av0]
acceptanceratesvsbeta = array(acceptanceratesvsbeta)

arstds = [std5, std2p5, std1, std0p75, std0p5, std0p4, std0p3, std0p2, std0p1, std0p05, std0p03, std0p01, std0]
arstds = array(arstds)

figure()
errorbar(betas, acceptanceratesvsbeta, yerr=arstds)
title('Acceptance rates for 10x10x10 fcc with $J=Dix=Diy=Diz=1$', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel('Acceptance rate', fontsize=20)
show()


    
    






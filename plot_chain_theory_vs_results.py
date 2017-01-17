from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys



# Numerical results
### L = 2
# J = 1, so betas = betas*J
# Open the file for reading

#infile = open("chain2_periodic_iso1_beta0to4_everybeta.txt", "r")
infile = open("test_everybeta.txt", "r")  # Just to compare methods

# Getting lists ready to store the data
betas     = []          # List of beta values
energies2 = []
energies2_stdv = []
energiessq2 = []
energiessq2_stdv = []
cvs2 = []
cvs2_stdv = []
mxs_av2 = []
mxs_stdv2 = []
mys_av2 = []
mys_stdv2 = []
mzs_av2 = []
mzs_stdv2 = []
acceptancerates2 = []
acceptancerates2_stdv = []
mxsqs_av2     = []
mxsqs_stdv2   = []
mysqs_av2     = []
mysqs_stdv2   = []
mzsqs_av2     = []
mzsqs_stdv2   = []
mxquads_av2   = []
mxquads_stdv2 = []
myquads_av2   = []
myquads_stdv2 = []
mzquads_av2   = []
mzquads_stdv2 = []


'''
beta " " energy_av " " E_stdv " " energy_sq_av " " Esq_stdv " " cv " " cv_stdv " " mx_av " "
mx_stdv " " my_av " " my_stdv " " mz_av " " mz_stdv " " ar_av " " ar_stdv;
'''

# Read the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        betas.append(float(words[0]))
        energies2.append(float(words[1]))
        energies2_stdv.append(float(words[2]))
        energiessq2.append(float(words[3]))
        energiessq2_stdv.append(float(words[4]))
        cvs2.append(float(words[5]))
        cvs2_stdv.append(float(words[6]))
        mxs_av2.append(float(words[7]))
        mxs_stdv2.append(float(words[8]))
        mys_av2.append(float(words[9]))
        mys_stdv2.append(float(words[10]))
        mzs_av2.append(float(words[11]))
        mzs_stdv2.append(float(words[12]))
        acceptancerates2.append(float(words[13]))
        acceptancerates2_stdv.append(float(words[14]))
        mxsqs_av2.append(float(words[15]))
        mxsqs_stdv2.append(float(words[16]))
        mysqs_av2.append(float(words[17]))
        mysqs_stdv2.append(float(words[18]))
        mzsqs_av2.append(float(words[19]))
        mzsqs_stdv2.append(float(words[20]))
        mxquads_av2.append(float(words[21]))
        mxquads_stdv2.append(float(words[22]))
        myquads_av2.append(float(words[23]))
        myquads_stdv2.append(float(words[24]))
        mzquads_av2.append(float(words[25]))
        mzquads_stdv2.append(float(words[26]))
        
# Remember to close the file
infile.close()

betas                 = array(betas)
energies2             = array(energies2)
energies2_stdv        = array(energies2_stdv)
energiessq2           = array(energiessq2)
energiessq2_stdv      = array(energiessq2_stdv)
cvs2                  = array(cvs2)#/2.            # Check this
cvs2_stdv             = array(cvs2_stdv)#/2.       # Check this
mxs_av2               = array(mxs_av2)
mxs_stdv2             = array(mxs_stdv2)
mys_av2               = array(mys_av2)
mys_stdv2             = array(mys_stdv2)
mzs_av2               = array(mzs_av2)
mzs_stdv2             = array(mzs_stdv2)
acceptancerates2      = array(acceptancerates2)
acceptancerates2_stdv = array(acceptancerates2_stdv)
mxsqs_av2     = array(mxsqs_av2)
mxsqs_stdv2   = array(mxsqs_stdv2)
mysqs_av2     = array(mysqs_av2)
mysqs_stdv2   = array(mysqs_stdv2)
mzsqs_av2     = array(mzsqs_av2)
mzsqs_stdv2   = array(mzsqs_stdv2)
mxquads_av2   = array(mxquads_av2)
mxquads_stdv2 = array(mxquads_stdv2)
myquads_av2   = array(myquads_av2)
myquads_stdv2 = array(myquads_stdv2)
mzquads_av2   = array(mzquads_av2)
mzquads_stdv2 = array(mzquads_stdv2)

# Theory. Since J =1, betas = betas*J 
'''
#Partition function. Perhaps I don't need this after all...
partition_function2 = zeros(len(betas))
partition_function2[0] = 16*pi**2 # Series expansion of sinh, inserting for beta*J = 0 afterwards
for i in range(1, len(betas)):
    partition_function2[i] = 8*pi**2/(betas[i])*sinh(2*betas[i])
'''
    
# Energy with J = 1
energy2_theory = 2*(1/(2*betas)-(cosh(2*betas)/sinh(2*betas)))
energy2_theory = zeros(len(betas))
for i in range(0,len(betas)):
    energy2_theory[i] =  2*( 1/(2*betas[i]) - (cosh(2*betas[i])/sinh(2*betas[i])))

# Heat capacity with J = 1
#k = 8.6173303e-5 # Boltzmann constant in eV K**(-1)
cvs2_theory = (1+4*(1-(1/sinh(2*betas)**2)))#*k

### Plotting
# Energy
figure()
errorbar(betas, energies2, yerr=energies2_stdv, label='Simulation')
hold('on')
plot(betas, energy2_theory, label='Theory')
title('Energy for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel('Energy', fontsize=20)
legend(loc='upper right')
show()

# Energy squared
figure()
errorbar(betas, energiessq2, yerr=energiessq2_stdv, label='Simulation')
title('Energy squared for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel('Energy', fontsize=20)
legend(loc='upper right')
show()

# Heat capacity
figure()
errorbar(betas, cvs2, yerr=cvs2_stdv, label='Simulation')
hold('on')
plot(betas, cvs2_theory)
title(r'Heat capacity $C_V$ for two particle chain', fontsize=16) # Specific?
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$C_V$', fontsize=20)
legend(loc='upper right')
show()

# Magnetization in the x-direction
figure()
errorbar(betas, mxs_av2, yerr=mxs_stdv2, label='Simulation')
title(r'Magnetization in the x-direction, $<m_x>$, for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$<m_x>$', fontsize=20)
legend(loc='upper right')
show()

# Magnetization in the y-direction
figure()
errorbar(betas, mys_av2, yerr=mys_stdv2, label='Simulation')
title(r'Magnetization in the y-direction, $<m_y>$, for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$<m_y>$', fontsize=20)
legend(loc='upper right')
show()

# Magnetization in the z-direction
figure()
errorbar(betas, mzs_av2, yerr=mzs_stdv2, label='Simulation')
title(r'Magnetization in the z-direction, $<m_z>$, for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$<m_z>$', fontsize=20)
legend(loc='upper right')
show()

# Magnetization in the x-direction
figure()
errorbar(betas, mxsqs_av2, yerr=mxsqs_stdv2, label='Simulation')
title(r'Magnetization squared in the x-direction, $<m^2_x>$, for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$<m^2_x>$', fontsize=20)
legend(loc='upper right')
show()

# Magnetization in the y-direction
figure()
errorbar(betas, mysqs_av2, yerr=mysqs_stdv2, label='Simulation')
title(r'Magnetization squared  in the y-direction, $<m^2_y>$, for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$<m^2_y>$', fontsize=20)
legend(loc='upper right')
show()

# Magnetization in the z-direction
figure()
errorbar(betas, mzsqs_av2, yerr=mzsqs_stdv2, label='Simulation')
title(r'Magnetization squared in the z-direction, $<m^2_z>$, for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$<m^2_z>$', fontsize=20)
legend(loc='upper right')
show()

# Magnetization in the x-direction
figure()
errorbar(betas, mxquads_av2, yerr=mxquads_stdv2, label='Simulation')
title(r'Quadruple magnetization in the x-direction, $<m^4_x>$, for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$<m^4_x>$', fontsize=20)
legend(loc='upper right')
show()

# Magnetization in the y-direction
figure()
errorbar(betas, myquads_av2, yerr=myquads_stdv2, label='Simulation')
title(r'Quadruple magnetization in the y-direction, $<m^4_y>$, for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$<m^4_y>$', fontsize=20)
legend(loc='upper right')
show()

# Magnetization in the z-direction
figure()
errorbar(betas, mzquads_av2, yerr=mzquads_stdv2, label='Simulation')
title(r'Quadruple magnetization  in the z-direction, $<m^4_z>$, for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$<m^4_z>$', fontsize=20)
legend(loc='upper right')
show()




# Acceptance rates
figure()
#errorbar(betas, acceptancerates2, yerr=acceptancerates2_stdv) #Weird things have happened
plot(betas, acceptancerates2)
title('Acceptance rates for two particle chain', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel('Acceptance rate', fontsize=20)
show()

'''
### Plotting
figure()
errorbar(betas, acceptancerates2, yerr=acceptancerates2_stdv, label='J=0.1')
hold('on')
errorbar(betas, acceptanceratesvsbeta0p5, yerr=arstds0p5, label='J=0.5')
hold('on')
errorbar(betas, acceptanceratesvsbeta1, yerr=arstds1, label='J=1')
hold('on')
errorbar(betas, acceptanceratesvsbeta5, yerr=arstds5, label='J=5')
hold('on')
errorbar(betas, acceptanceratesvsbeta10, yerr=arstds10, label='J=10')
hold('on')
errorbar(betas, acceptanceratesvsbeta50, yerr=arstds50, label='J=50')
hold('on')
title('Acceptance rates for isotropic 10x10x10 fcc', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel('Acceptance rate', fontsize=20)
legend(loc='upper right')
show()

'''



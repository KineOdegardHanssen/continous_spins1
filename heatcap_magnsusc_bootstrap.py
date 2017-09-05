from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys
from random import randint

#############################################################################
# Extracting heat capacity as a function of temperature, etc.

def get_Bootstrapresults(filename, L):
    N = L*L*L
    infile = open(filename, "r")    
    # Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
    
    # Sorting arrays
    betasA                 = []          # List of beta values
    binofbeta_firstindexA  = []          # I guess nested lists are a bit of a hassle
    binofbeta_lastindexA   = []          # Could also just feed the number of bins in. But, the less such input, the better, probably
    # Is the above unsafe? Should I just find the no of bins and add instead, a la, i*nbins+something?
    
    # Other quantities
    energy_bavsA      = []          # Index 1
    energy_sq_bavsA   = []          # Index 2
    hc_bavsA          = []          # Index 3
    mxsc_bavsA        = []          # Index 4
    mxsc_abs_bavsA    = []          # Index 5
    mxssqc_bavsA      = []          # Index 6
    mxsquadc_bavsA    = []          # Index 7
    mysc_bavsA        = []          # Index 8
    mysc_abs_bavsA    = []          # Index 9
    myssqc_bavsA      = []          # Index 10
    mysquadc_bavsA    = []          # Index 11
    mzsc_bavsA        = []          # Index 12
    mzsc_abs_bavsA    = []          # Index 13
    mzssqc_bavsA      = []          # Index 14
    mzsquadc_bavsA    = []          # Index 15
    
    # Don't think I need the magnetizations of q=(0,0,0), but we'll see...
    # Read the rest of the lines
    lines = infile.readlines()  # This works, I have checked it
    
    i = 0 # Some sort of counter to get the indices right
    betabefore = 0 # Test the betas so we can gather the bins.
    binofbeta_firstindexA.append(0) 
    # Should I have some list over the indices of the betas belonging to each bin?
    
    # Getting data from the file
    for line in lines:
        words = line.split()
        if len(words) != 0:
            # Betas
            beta = float(words[0])
            if beta != betabefore:
                betasA.append(beta)
                betabefore = beta
                if i != 0:
                    binofbeta_firstindexA.append(i)
                    binofbeta_lastindexA.append(i-1)
            # energy
            en = float(words[1])
            energy_bavsA.append(en)          # Index 1
            # energy squared
            en = float(words[2])
            energy_sq_bavsA.append(en)       # Index 2
            # heat capacity
            en = float(words[3])
            hc_bavsA.append(en)              # Index 3
            # <m_x(q)>
            en = float(words[4])
            mxsc_bavsA.append(en)            # Index 4
            # <|m_x(q)|>
            en = float(words[5])
            mxsc_abs_bavsA.append(en)        # Index 5
            # <m^2_x(q)>
            en = float(words[6])
            mxssqc_bavsA.append(en)          # Index 6
            # <m^4_x(q)>
            en = float(words[7])
            mxsquadc_bavsA.append(en)        # Index 7
            # <m_y(q)>
            en = float(words[8])
            mysc_bavsA.append(en)            # Index 8
            # <|m_y(q)|>
            en = float(words[9])
            mysc_abs_bavsA.append(en)        # Index 9
            # <m^2_y(q)>
            en = float(words[10])
            myssqc_bavsA.append(en)          # Index 10
            # <m^4_y(q)>
            en = float(words[11])
            mysquadc_bavsA.append(en)        # Index 11
            # <m_z(q)>
            en = float(words[12])
            mzsc_bavsA.append(en)            # Index 12
            # <|m_z(q)|>
            en = float(words[13])
            mzsc_abs_bavsA.append(en)        # Index 13
            # <m^2_z(q)>
            en = float(words[14])
            mzssqc_bavsA.append(en)          # Index 14
            # <m^4_z(q)>
            en = float(words[15])
            mzsquadc_bavsA.append(en)        # Index 15
        i += 1 # Increasing the counter
        
    binofbeta_lastindexA.append(i-1)  # The index has been updated one time too many
    
    # Remember to close the file
    infile.close()
    
    no_of_bins_each_betaA = binofbeta_lastindexA[0]+1 # Due to the indexing in Python + we don't allow
    no_of_betasA          = len(betasA)               # The number of different temperatures ### Should I convert to T now?
    no_of_Bootstrapruns   = no_of_bins_each_betaA     # Probably OK this way... Good to be flexible...
    nbins                 = no_of_bins_each_betaA
    
    betasA = array(betasA) # Not sure if this is neccessary, but maybe...
    # Magnetic susceptibility per spin:
    magnsuscx_vec = zeros(no_of_betasA); magnsuscy_vec = zeros(no_of_betasA); magnsuscz_vec = zeros(no_of_betasA);# Not quite sure how to order this...
    magnsuscx_stdv_vec = zeros(no_of_betasA); magnsuscy_stdv_vec = zeros(no_of_betasA); magnsuscz_stdv_vec = zeros(no_of_betasA);# Not quite sure how to order this...    
    
    magnsuscx_br = zeros(no_of_Bootstrapruns)
    magnsuscy_br = zeros(no_of_Bootstrapruns)
    magnsuscz_br = zeros(no_of_Bootstrapruns)
    
    # Total heat capacity:
    heatcap_vec = zeros(no_of_betasA)
    heatcap_stdv_vec = zeros(no_of_betasA)
    
    heatcap_br = zeros(no_of_Bootstrapruns)
    
    ### Oy, vey, all the loops... # Need one for several bootstrap runs too.
    ### for different numbers of bins for different beta
    for i in range(0, no_of_betasA): # For doing it for every temperature
        for k in range(0, no_of_Bootstrapruns): # Want to repeat the random selection of bins a number of times
            ### Reset all quantities I want to find by bootstrap        
            mx = 0; mxabs = 0; mx2 = 0; mx4 = 0; my = 0; myabs = 0; my2 = 0; my4 = 0; mz = 0; mzabs = 0; mz2 = 0; mz4 = 0;
            eav = 0; esqav = 0; hcav = 0;
            for j in range(0, nbins): # For finding the averages
                ###Draw a random integer in [binofbeta_firstindexA[i], binofbeta_lastindexA[i]].
                n = randint(binofbeta_firstindexA[i], binofbeta_lastindexA[i])  # Draw a random number n times
                ###Extract O(n), add to average function
                eav += energy_bavsA[n]; esqav += energy_sq_bavsA[n]; hcav += hc_bavsA[n]
                mx += mxsc_bavsA[n]; mxabs += mxsc_abs_bavsA[n]; mx2 += mxssqc_bavsA[n]; mx4 += mxsquadc_bavsA[n];
                my += mysc_bavsA[n]; myabs += mysc_abs_bavsA[n]; my2 += myssqc_bavsA[n]; my4 += mysquadc_bavsA[n];
                mz += mzsc_bavsA[n]; mzabs += mzsc_abs_bavsA[n]; mz2 += mzssqc_bavsA[n]; mz4 += mzsquadc_bavsA[n];
            ### Divide the average of O(n) by no_of_bins_each_beta to ACTUALLY make the average
            eav = eav/nbins; esqav = esqav/nbins; hcav = hcav/nbins
            mx = mx/nbins; mxabs = mxabs/nbins; mx2 = mx2/nbins; mx4 = mx4/nbins;
            my = my/nbins; myabs = myabs/nbins; my2 = my2/nbins; my4 = my4/nbins;
            mz = mz/nbins; mzabs = mzabs/nbins; mz2 = mz2/nbins; mz4 = mz4/nbins;
            ### Probably store it in a list. Ya know, for standard deviations and such
            ### Magnetic susceptibility per spin
            magnsuscx_this = betasA[i]*(mx2-(mxabs*mxabs))*N # OFS said I should use the absolute magn to the first power. AF, prob 
            magnsuscy_this = betasA[i]*(my2-(myabs*myabs))*N
            magnsuscz_this = betasA[i]*(mz2-(mzabs*mzabs))*N # But what should I do with them
            '''
            magnsuscx_this = betasA[i]*(mx2-(mx*mx)) # OFS said I should use the absolute magn to the first power. AF, prob 
            magnsuscy_this = betasA[i]*(my2-(my*my))
            magnsuscz_this = betasA[i]*(mz2-(mz*mz)) # But what should I do with them
            '''
            # Finding the average magnetic susceptibility
            magnsuscx_vec[i] += magnsuscx_this
            magnsuscy_vec[i] += magnsuscy_this
            magnsuscz_vec[i] += magnsuscz_this
            # Storing for every Bootstrap-iteration
            magnsuscx_br[k] = magnsuscx_this
            magnsuscy_br[k] = magnsuscy_this
            magnsuscz_br[k] = magnsuscz_this
            ### Total heat capacity
            heatcap_this = betasA[i]**2*(esqav-(eav*eav))/N
            # Finding the average heat capacity
            heatcap_vec[i] += heatcap_this
            # Storing for every Bootstrap-iteration
            heatcap_br[k] = heatcap_this
            ### Find other stuff as well ??
            ### Do a fitting of the results. To a line, a quadr. func., a pol. of the third degree. Check this. OFS suggested quadr
        # Treating the magnetic susceptibility data
        # Finding the average
        magnsuscx_vec[i] = magnsuscx_vec[i]/no_of_Bootstrapruns
        magnsuscy_vec[i] = magnsuscy_vec[i]/no_of_Bootstrapruns
        magnsuscz_vec[i] = magnsuscz_vec[i]/no_of_Bootstrapruns  
        # Finding the standard deviation:
        for k in range(0, no_of_Bootstrapruns):
            magnsuscx_stdv_vec[i] += (magnsuscx_br[k]-magnsuscx_vec[i])*(magnsuscx_br[k]-magnsuscx_vec[i])
            magnsuscy_stdv_vec[i] += (magnsuscy_br[k]-magnsuscy_vec[i])*(magnsuscy_br[k]-magnsuscy_vec[i])
            magnsuscz_stdv_vec[i] += (magnsuscz_br[k]-magnsuscz_vec[i])*(magnsuscz_br[k]-magnsuscz_vec[i])
        magnsuscx_stdv_vec[i] = sqrt(magnsuscx_stdv_vec[i]/(no_of_Bootstrapruns*(no_of_Bootstrapruns-1)))
        magnsuscy_stdv_vec[i] = sqrt(magnsuscy_stdv_vec[i]/(no_of_Bootstrapruns*(no_of_Bootstrapruns-1)))
        magnsuscz_stdv_vec[i] = sqrt(magnsuscz_stdv_vec[i]/(no_of_Bootstrapruns*(no_of_Bootstrapruns-1)))
        
        # Treating the heat capacity data
        # Finding the average
        heatcap_vec[i] = heatcap_vec[i]/no_of_Bootstrapruns
        # Finding the standard deviation
        for k in range(0, no_of_Bootstrapruns): # Newman & Barkema gave another expression for the error in the Bootstrap chapter.
            heatcap_stdv_vec[i] += (heatcap_br[k]-heatcap_vec[i])*(heatcap_br[k]-heatcap_vec[i])
        heatcap_stdv_vec[i] = sqrt(heatcap_stdv_vec[i]/(no_of_Bootstrapruns*(no_of_Bootstrapruns-1)))
    return betasA, magnsuscx_vec, magnsuscy_vec, magnsuscz_vec,magnsuscx_stdv_vec, magnsuscy_stdv_vec, magnsuscz_stdv_vec, heatcap_vec, heatcap_stdv_vec


#################### This is where we give our input #######################
LA = 6
LB = 8
LC = 10
LD = 12
LE = 14

# Jensen's parameters  ## Please keep and comment out when using other param.s. Saving a lot of bother...
filenameA = "fcc6x6x6yopen_beta0p775to0p779_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameB = "fcc8x8x8yopen_beta0p775to0p779_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameC = "fcc10x10x10yopen_beta0p775to0p779_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameD = "fcc12x12x12yopen_beta0p775to0p779_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameE = "fcc14x14x14yopen_beta0p775to0p779_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
##################                                   #######################

# Starting the Bootstrap process
betaarrayA, magnsuscx_vecA, magnsuscy_vecA, magnsuscz_vecA, magnsuscx_stdv_vecA, magnsuscy_stdv_vecA, magnsuscz_stdv_vecA, heatcap_vecA, heatcap_stdv_vecA = get_Bootstrapresults(filenameA, LA)
betaarrayB, magnsuscx_vecB, magnsuscy_vecB, magnsuscz_vecB, magnsuscx_stdv_vecB, magnsuscy_stdv_vecB, magnsuscz_stdv_vecB, heatcap_vecB, heatcap_stdv_vecB = get_Bootstrapresults(filenameB, LB)
betaarrayC, magnsuscx_vecC, magnsuscy_vecC, magnsuscz_vecC, magnsuscx_stdv_vecC, magnsuscy_stdv_vecC, magnsuscz_stdv_vecC, heatcap_vecC, heatcap_stdv_vecC = get_Bootstrapresults(filenameC, LC)
betaarrayD, magnsuscx_vecD, magnsuscy_vecD, magnsuscz_vecD, magnsuscx_stdv_vecD, magnsuscy_stdv_vecD, magnsuscz_stdv_vecD, heatcap_vecD, heatcap_stdv_vecD = get_Bootstrapresults(filenameD, LD)
betaarrayE, magnsuscx_vecE, magnsuscy_vecE, magnsuscz_vecE, magnsuscx_stdv_vecE, magnsuscy_stdv_vecE, magnsuscz_stdv_vecE, heatcap_vecE, heatcap_stdv_vecE = get_Bootstrapresults(filenameE, LE)

'''
# "Normalizing" the functions (according to page 17 in Barkema Newman)
magnsuscx_vecA = magnsuscx_vecA*NA*NA; magnsuscy_vecA = magnsuscy_vecA*NA*NA; magnsuscz_vecA = magnsuscz_vecA*NA
magnsuscx_stdv_vecA = magnsuscx_stdv_vecA*NA*NA; magnsuscy_stdv_vecA = magnsuscy_stdv_vecA*NA*NA; magnsuscz_stdv_vecA = magnsuscz_stdv_vecA*NA*NA

magnsuscx_vecB = magnsuscx_vecB*NB*NB; magnsuscy_vecB = magnsuscy_vecB*NB*NB; magnsuscz_vecB = magnsuscz_vecB*NB*NB
magnsuscx_stdv_vecB = magnsuscx_stdv_vecB*NB*NB; magnsuscy_stdv_vecB = magnsuscy_stdv_vecB*NB*NB; magnsuscz_stdv_vecB = magnsuscz_stdv_vecB*NB*NB

magnsuscx_vecC = magnsuscx_vecC*NC; magnsuscy_vecC = magnsuscy_vecC*NC*NC; magnsuscz_vecC = magnsuscz_vecC*NC*NC
magnsuscx_stdv_vecC = magnsuscx_stdv_vecC*NC*NC; magnsuscy_stdv_vecC = magnsuscy_stdv_vecC*NC*NC; magnsuscz_stdv_vecC = magnsuscz_stdv_vecC*NC*NC

magnsuscx_vecE = magnsuscx_vecE*NE*NE; magnsuscy_vecE = magnsuscy_vecE*NE*NE; magnsuscz_vecA = magnsuscz_vecE*NE*NE
magnsuscx_stdv_vecE = magnsuscx_stdv_vecE*NE*NE; magnsuscy_stdv_vecE = magnsuscy_stdv_vecE*NE*NE; magnsuscz_stdv_vecE = magnsuscz_stdv_vecE*NE*NE

magnsuscx_vecD = magnsuscx_vecD*ND*ND; magnsuscy_vecD = magnsuscy_vecD*ND*ND; magnsuscz_vecD = magnsuscz_vecD*ND*ND
magnsuscx_stdv_vecD = magnsuscx_stdv_vecD*ND*ND; magnsuscy_stdv_vecD = magnsuscy_stdv_vecD*ND*ND; magnsuscz_stdv_vecD = magnsuscz_stdv_vecD*ND*ND

heatcap_vecA = heatcap_vecA/float(NA); heatcap_stdv_vecA = heatcap_stdv_vecA/float(NA)
heatcap_vecB = heatcap_vecA/float(NB); heatcap_stdv_vecB = heatcap_stdv_vecB/float(NB)
heatcap_vecC = heatcap_vecA/float(NC); heatcap_stdv_vecC = heatcap_stdv_vecC/float(NC)
heatcap_vecD = heatcap_vecA/float(ND); heatcap_stdv_vecD = heatcap_stdv_vecD/float(ND)
heatcap_vecE = heatcap_vecA/float(NE); heatcap_stdv_vecE = heatcap_stdv_vecE/float(NE)
'''

# Plotting
figure(figsize=(6,5))
errorbar(betaarrayA, magnsuscz_vecA, yerr=magnsuscz_stdv_vecA, capsize=2, label='L=%i'%LA)
hold('on')
errorbar(betaarrayB, magnsuscz_vecB, yerr=magnsuscz_stdv_vecB, capsize=2, label='L=%i'%LB)
errorbar(betaarrayC, magnsuscz_vecC, yerr=magnsuscz_stdv_vecC, capsize=2, label='L=%i'%LC)
errorbar(betaarrayD, magnsuscz_vecD, yerr=magnsuscz_stdv_vecD, capsize=2, label='L=%i'%LD)
errorbar(betaarrayE, magnsuscz_vecE, yerr=magnsuscz_stdv_vecE, capsize=2, label='L=%i'%LE)
title(r'Magnetic susceptibility vs $\beta$ for different L', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$\chi_M/N$', fontsize=20)
legend(loc="upper left")
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
show()

figure(figsize=(6,5))
errorbar(betaarrayA, heatcap_vecA, yerr=heatcap_stdv_vecA, capsize=2, label='L=%i'%LA)
hold('on')
errorbar(betaarrayB, heatcap_vecB, yerr=heatcap_stdv_vecB, capsize=2, label='L=%i'%LB)
errorbar(betaarrayC, heatcap_vecC, yerr=heatcap_stdv_vecC, capsize=2, label='L=%i'%LC)
errorbar(betaarrayD, heatcap_vecD, yerr=heatcap_stdv_vecD, capsize=2, label='L=%i'%LD)
errorbar(betaarrayE, heatcap_vecE, yerr=heatcap_stdv_vecE, capsize=2, label='L=%i'%LE)
title(r'Heat capacity vs $\beta$ for different L', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$c_V$', fontsize=20)
legend(loc="upper left")
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
show()

'''
# Magn susc in all dir

figure(figsize=(6,5))
errorbar(betaarrayA, magnsuscx_vecA+magnsuscy_vecA+magnsuscz_vecA, yerr=magnsuscx_stdv_vecA+magnsuscy_stdv_vecA+magnsuscz_stdv_vecA, capsize=2, label='L=%i'%LA)
hold('on')
errorbar(betaarrayB, magnsuscx_vecB+magnsuscy_vecB+magnsuscz_vecB, yerr=magnsuscx_stdv_vecB+magnsuscy_stdv_vecB+magnsuscz_stdv_vecB, capsize=2, label='L=%i'%LB)
errorbar(betaarrayC, magnsuscx_vecC+magnsuscy_vecC+magnsuscz_vecC, yerr=magnsuscx_stdv_vecC+magnsuscy_stdv_vecC+magnsuscz_stdv_vecC, capsize=2, label='L=%i'%LC)
errorbar(betaarrayD, magnsuscx_vecD+magnsuscy_vecD+magnsuscz_vecD, yerr=magnsuscx_stdv_vecD+magnsuscy_stdv_vecD+magnsuscz_stdv_vecD, capsize=2, label='L=%i'%LD)
errorbar(betaarrayE, magnsuscx_vecE+magnsuscy_vecE+magnsuscz_vecE, yerr=magnsuscx_stdv_vecE+magnsuscy_stdv_vecE+magnsuscz_stdv_vecE, capsize=2, label='L=%i'%LE)
title(r'Magnetic susceptibility vs $\beta$ for different L', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$\chi_M$', fontsize=20)
legend(loc="upper left")
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
show()
'''

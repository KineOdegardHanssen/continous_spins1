from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys
from random import randint

'''
# Test of rng
for i in range(8):
    print randint(0,3)
'''

# I will extract only info about Tc from the Binder cumulant first

def givequadratic(a, b, c, x): # Should this really be quadratic? I look at a small interval...
    return a*x**2+b*x+c

def findintersection(beta, f1, f2):
    diff = f1 - f2
    for i in range(len(diff)):
        if (i+1)<len(diff): # Just a safeguard. We're doomed if this happens anyway
            if (diff[i] == 0) or (diff[i]*diff[i+1] < 0):
                return beta[i]
    print "NB! Intersection not found!"
    return beta[0]

def extract_crossing(filenameA, filenameB): # Extract the beta value of the crossing between two graphs
    infileA = open(filenameA, "r")
    # Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
    
    # Sorting arrays
    betasA                 = []          # List of beta values
    binofbeta_firstindexA  = []          # I guess nested lists are a bit of a hassle
    binofbeta_lastindexA   = []          # Could also just feed the number of bins in. But, the less such input, the better, probably
    # Is the above unsafe? Should I just find the no of bins and add instead, a la, i*nbins+something?
    
    # Other quantities
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
    lines = infileA.readlines()  # This works, I have checked it
    
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
    infileA.close()
    
    no_of_bins_each_betaA = binofbeta_lastindexA[0]+1 # Due to the indexing in Python + we don't allow
    no_of_betasA          = len(betasA)               # The number of different temperatures ### Should I convert to T now?
    no_of_BootstraprunsA  = no_of_bins_each_betaA     # Probably OK this way... Good to be flexible...
    nbinsA                = no_of_bins_each_betaA
    
    betasA = array(betasA) # Not sure if this is neccessary, but maybe...
    
    # Then B
    
    infileB = open(filenameB, "r")
    # Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
    
    # Sorting arrays
    betasB                 = []          # List of beta values
    binofbeta_firstindexB  = []          # I guess nested lists are a bit of a hassle
    binofbeta_lastindexB   = []          # Could also just feed the number of bins in. But, the less such input, the better, probably
    # Is the above unsafe? Should I just find the no of bins and add instead, a la, i*nbins+something?
    
    # Other quantities
    mxsc_bavsB        = []          # Index 4
    mxsc_abs_bavsB    = []          # Index 5
    mxssqc_bavsB      = []          # Index 6
    mxsquadc_bavsB    = []          # Index 7
    mysc_bavsB        = []          # Index 8
    mysc_abs_bavsB    = []          # Index 9
    myssqc_bavsB      = []          # Index 10
    mysquadc_bavsB    = []          # Index 11
    mzsc_bavsB        = []          # Index 12
    mzsc_abs_bavsB    = []          # Index 13
    mzssqc_bavsB      = []          # Index 14
    mzsquadc_bavsB    = []          # Index 15
    
    # Don't think I need the magnetizations of q=(0,0,0), but we'll see...
    # Read the rest of the lines
    lines = infileB.readlines()  # This works, I have checked it
    
    i = 0 # Some sort of counter to get the indices right
    betabefore = 0 # Test the betas so we can gather the bins.
    binofbeta_firstindexB.append(0) 
    # Should I have some list over the indices of the betas belonging to each bin?
    
    # Getting data from the file
    for line in lines:
        words = line.split()
        if len(words) != 0:
            # Betas
            beta = float(words[0])
            if beta != betabefore:
                betasB.append(beta)
                betabefore = beta
                if i != 0:
                    binofbeta_firstindexB.append(i)
                    binofbeta_lastindexB.append(i-1)
            # <m_x(q)>
            en = float(words[4])
            mxsc_bavsB.append(en)            # Index 4
            # <|m_x(q)|>
            en = float(words[5])
            mxsc_abs_bavsB.append(en)        # Index 5
            # <m^2_x(q)>
            en = float(words[6])
            mxssqc_bavsB.append(en)          # Index 6
            # <m^4_x(q)>
            en = float(words[7])
            mxsquadc_bavsB.append(en)        # Index 7
            # <m_y(q)>
            en = float(words[8])
            mysc_bavsB.append(en)            # Index 8
            # <|m_y(q)|>
            en = float(words[9])
            mysc_abs_bavsB.append(en)        # Index 9
            # <m^2_y(q)>
            en = float(words[10])
            myssqc_bavsB.append(en)          # Index 10
            # <m^4_y(q)>
            en = float(words[11])
            mysquadc_bavsB.append(en)        # Index 11
            # <m_z(q)>
            en = float(words[12])
            mzsc_bavsB.append(en)            # Index 12
            # <|m_z(q)|>
            en = float(words[13])
            mzsc_abs_bavsB.append(en)        # Index 13
            # <m^2_z(q)>
            en = float(words[14])
            mzssqc_bavsB.append(en)          # Index 14
            # <m^4_z(q)>
            en = float(words[15])
            mzsquadc_bavsB.append(en)        # Index 15
        i += 1 # Increasing the counter
        
    binofbeta_lastindexB.append(i-1)  # The index has been updated one time too many
    
    # Remember to close the file
    infileB.close()
    
    no_of_bins_each_betaB = binofbeta_lastindexB[0]+1 # Due to the indexing in Python + we don't allow
    no_of_betasB          = len(betasB)               # The number of different temperatures ### Should I convert to T now?
    no_of_BootstraprunsB  = no_of_bins_each_betaB     # Probably OK this way... Good to be flexible...
    nbinsB                = no_of_bins_each_betaB
    
    betasA = array(betasA) # Not sure if this is neccessary, but maybe...
    
    intersectionbeta_av = 0
    intersectionbetas = zeros(no_of_BootstraprunsB)
    for k in range(0, no_of_BootstraprunsA): # Want to repeat the random selection of bins a number of times
        bindervalueszcompA = []
        bindervalueszcompB = []
        for i in range(0, no_of_betasA): # For doing it for every temperature
            ### Reset all quantities I want to find by bootstrap        
            # From A
            mxA = 0; mxabsA = 0; mx2A = 0; mx4A = 0;
            myA = 0; myabsA = 0; my2A = 0; my4A = 0;
            mzA = 0; mzabsA = 0; mz2A = 0; mz4A = 0;
            # From B
            mxB = 0; mxabsB = 0; mx2B = 0; mx4B = 0;
            myB = 0; myabsB = 0; my2B = 0; my4B = 0;
            mzB = 0; mzabsB = 0; mz2B = 0; mz4B = 0;
            for j in range(0, nbinsA): # For finding the averages
                ###Draw a random integer in [binofbeta_firstindexA[i], binofbeta_lastindexA[i]].
                n = randint(binofbeta_firstindexA[i], binofbeta_lastindexA[i])  # Draw a random number n times
                ###Extract O(n), add to average function
                # From A
                mxA += mxsc_bavsA[n]; mxabsA += mxsc_abs_bavsA[n]; mx2A += mxssqc_bavsA[n]; mx4A += mxsquadc_bavsA[n];
                myA += mysc_bavsA[n]; myabsA += mysc_abs_bavsA[n]; my2A += myssqc_bavsA[n]; my4A += mysquadc_bavsA[n];
                mzA += mzsc_bavsA[n]; mzabsA += mzsc_abs_bavsA[n]; mz2A += mzssqc_bavsA[n]; mz4A += mzsquadc_bavsA[n];
                # From B
                mxB += mxsc_bavsB[n]; mxabsB += mxsc_abs_bavsB[n]; mx2B += mxssqc_bavsB[n]; mx4B += mxsquadc_bavsB[n];
                myB += mysc_bavsB[n]; myabsB += mysc_abs_bavsB[n]; my2B += myssqc_bavsB[n]; my4B += mysquadc_bavsB[n];
                mzB += mzsc_bavsB[n]; mzabsB += mzsc_abs_bavsB[n]; mz2B += mzssqc_bavsB[n]; mz4B += mzsquadc_bavsB[n];
            ### Divide the average of O(n) by no_of_bins_each_beta to ACTUALLY make the average
            # From A
            mxA = mxA/nbinsA; mxabsA = mxabsA/nbinsA; mx2A = mx2A/nbinsA; mx4A = mx4A/nbinsA;
            myA = myA/nbinsA; myabsA = myabsA/nbinsA; my2A = my2A/nbinsA; my4A = my4A/nbinsA;
            mzA = mzA/nbinsA; mzabsA = mzabsA/nbinsA; mz2A = mz2A/nbinsA; mz4A = mz4A/nbinsA;
            # From B
            mxB = mxB/nbinsA; mxabsB = mxabsB/nbinsA; mx2B = mx2B/nbinsA; mx4B = mx4B/nbinsA;
            myB = myB/nbinsA; myabsB = myabsB/nbinsA; my2B = my2B/nbinsA; my4B = my4B/nbinsA;
            mzB = mzB/nbinsA; mzabsB = mzabsB/nbinsA; mz2B = mz2B/nbinsA; mz4B = mz4B/nbinsA;
            ### Probably store it in a list. Ya know, for standard deviations and such
            ### Find the Binder cumulant
            # For A
            thisbinderzA = 1 - (mz4A/(3.*mz2A*mz2A))
            bindervalueszcompA.append(thisbinderzA)
            # For B
            thisbinderzB = 1 - (mz4B/(3.*mz2B*mz2B))
            bindervalueszcompB.append(thisbinderzB)
            
            ### Do a fitting of the results. To a line, a quadr. func., a pol. of the third degree. Check this. OFS suggested quadr
        # For A
        bindervalueszcompA = array(bindervalueszcompA) 
        fitvectorzcompA = polyfit(betasA, bindervalueszcompA, 2) # Fits the function points to a quadratic polynomial
        azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]
        # For B
        bindervalueszcompB = array(bindervalueszcompB) 
        fitvectorzcompB = polyfit(betasB, bindervalueszcompB, 2) # Fits the function points to a quadratic polynomial
        azB = fitvectorzcompB[0]; bzB = fitvectorzcompB[1]; czB = fitvectorzcompB[2]
        
        nfbetas = 1000
        fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
        
        fA = givequadratic(azA, bzA, czA, fbetas) # Should this really be quadratic? I look at a small interval...
        fB = givequadratic(azB, bzB, czB, fbetas)
        crossingbeta = findintersection(fbetas, fA, fB)
        intersectionbeta_av += crossingbeta
        intersectionbetas[k] = crossingbeta
        
    intersectionbeta_av = intersectionbeta_av/no_of_BootstraprunsA
    intersectionbeta_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        intersectionbeta_stddev += (intersectionbeta_av-intersectionbetas[i])*(intersectionbeta_av-intersectionbetas[i])
    intersectionbeta_stddev = sqrt(intersectionbeta_stddev/(no_of_BootstraprunsA*(no_of_BootstraprunsA-1)))
    return intersectionbeta_av, intersectionbeta_stddev


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

ib_avA, ib_stdvA = extract_crossing(filenameA, filenameB)
ib_avB, ib_stdvB = extract_crossing(filenameB, filenameC)
ib_avC, ib_stdvC = extract_crossing(filenameC, filenameD)
ib_avD, ib_stdvD = extract_crossing(filenameD, filenameE)

crossingbetas      = [ib_avA, ib_avB, ib_avC, ib_avD]
crossingbetas_stdv = [ib_stdvA, ib_stdvB, ib_stdvC, ib_stdvD]
vec1dL             = [1./LA, 1./LB, 1./LC, 1./LD]

crossingbetas      = array(crossingbetas)
crossingbetas_stdv = array(crossingbetas_stdv)
vec1dL             = array(vec1dL)

print "Min crossingbetas", min(crossingbetas)
print "Max crossingbetas", max(crossingbetas)

print "Min crossingbetas*0.98", min(crossingbetas)*0.999
print "Max crossingbetas*1.02", max(crossingbetas)*1.001

figure(figsize=(6,5))
errorbar(vec1dL, crossingbetas, yerr=crossingbetas_stdv, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$1/L$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, (max(vec1dL)*1.02), min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()




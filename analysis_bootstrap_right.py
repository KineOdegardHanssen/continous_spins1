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

def giveline(a,b,x):
    return a*x+b

def givequadratic(a, b, c, x): # Should this really be quadratic? I look at a small interval...
    return a*x**2+b*x+c
    
def givecubic(a,b,c,d,x):
    return a*x**3+b*x**2+c*x+d

def givequadruple(a,b,c,d,e,x):
    return a*x**4+b*x**3+c*x**2+d*x+e
    
def givepowerlaw(k,n,x):
    return k*x**n

def findintersection(beta, f1, f2):
    diff = f1 - f2
    for i in range(len(diff)):
        if (i+1)<len(diff): # Just a safeguard. We're doomed if this happens anyway
            if (diff[i] == 0) or (diff[i]*diff[i+1] < 0):
                return  beta[i], f1[i]
    print "NB! Intersection not found!"
    return beta[0], f1[0]

def extract_crossing(filenameA, filenameB, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL): # Extract the beta value of the crossing between two graphs
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
    
    nbinsA = zeros(len(betasA))
    for i in range(0,len(betasA)):
        nbinsA[i] = binofbeta_lastindexA[i]-binofbeta_firstindexA[i]+1
    
    no_of_bins_each_betaA = binofbeta_lastindexA[0]+1 # Due to the indexing in Python + we don't allow
    no_of_betasA          = len(betasA)               # The number of different temperatures ### Should I convert to T now?
    no_of_BootstraprunsA  = no_of_bins_each_betaA     # Probably OK this way... Good to be flexible...
    
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
    if toplotornot==0:
        plotter = 1
    else:
        plotter = 0
    
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
    
    nbinsB = zeros(len(betasB))
    for i in range(0,len(betasB)):
        nbinsB[i] = binofbeta_lastindexB[i]-binofbeta_firstindexB[i]+1
    
    no_of_bins_each_betaB = binofbeta_lastindexB[0]+1 # Due to the indexing in Python + we don't allow
    no_of_betasB          = len(betasB)               # The number of different temperatures ### Should I convert to T now?
    no_of_BootstraprunsB  = no_of_bins_each_betaB     # Probably OK this way... Good to be flexible...
    
    if cutit==1:   # If we have simulations for a larger temperature interval than we wish
        betamin = cutlower
        betamax = cutupper
        # For A
        smallbetaarray = []
        binofbeta_firstindex_thosewewant = []
        binofbeta_lastindex_thosewewant = []
        nbins_wewant = []
        no_of_betas_wewant= 0
        for i in range(0, no_of_betasA):
            if (betasA[i]<betamax and betasA[i]>betamin):
                smallbetaarray.append(betasA[i])
                nbins_wewant.append(nbinsA[i])
                binofbeta_firstindex_thosewewant.append(binofbeta_firstindexA[i])
                binofbeta_lastindex_thosewewant.append(binofbeta_lastindexA[i])
                no_of_betas_wewant+=1
        no_of_betasA = no_of_betas_wewant
        betasA = smallbetaarray
        nbinsA = nbins_wewant
        binofbeta_firstindexA = binofbeta_firstindex_thosewewant
        binofbeta_lastindexA = binofbeta_lastindex_thosewewant
        # For B # Well, A and B must be of the same size, I think, but this can't hurt
        smallbetaarray = []
        binofbeta_firstindex_thosewewant = []
        binofbeta_lastindex_thosewewant = []
        no_of_betas_wewant= 0
        for i in range(0, no_of_betasB):
            if (betasB[i]<betamax and betasB[i]>betamin):
                smallbetaarray.append(betasB[i])
                nbins_wewant.append(nbinsB[i])
                binofbeta_firstindex_thosewewant.append(binofbeta_firstindexB[i])
                binofbeta_lastindex_thosewewant.append(binofbeta_lastindexB[i])
                no_of_betas_wewant+=1
        no_of_betasB = no_of_betas_wewant
        betasB = smallbetaarray
        nbinsB = nbins_wewant
        binofbeta_firstindexB = binofbeta_firstindex_thosewewant
        binofbeta_lastindexB = binofbeta_lastindex_thosewewant
    betasA = array(betasA) # Not sure if this is neccessary, but maybe...
    
    intersectionbeta_av = 0
    intersectionbetas = zeros(no_of_BootstraprunsB)
    intersectionU_av = 0
    intersectionUs = zeros(no_of_BootstraprunsB)
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
            nbinsthis = int(nbinsA[i])
            for j in range(0, nbinsthis): # For finding the averages
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
            mxA = mxA/nbinsthis; mxabsA = mxabsA/nbinsthis; mx2A = mx2A/nbinsthis; mx4A = mx4A/nbinsthis;
            myA = myA/nbinsthis; myabsA = myabsA/nbinsthis; my2A = my2A/nbinsthis; my4A = my4A/nbinsthis;
            mzA = mzA/nbinsthis; mzabsA = mzabsA/nbinsthis; mz2A = mz2A/nbinsthis; mz4A = mz4A/nbinsthis;
            # From B
            mxB = mxB/nbinsthis; mxabsB = mxabsB/nbinsthis; mx2B = mx2B/nbinsthis; mx4B = mx4B/nbinsthis;
            myB = myB/nbinsthis; myabsB = myabsB/nbinsthis; my2B = my2B/nbinsthis; my4B = my4B/nbinsthis;
            mzB = mzB/nbinsthis; mzabsB = mzabsB/nbinsthis; mz2B = mz2B/nbinsthis; mz4B = mz4B/nbinsthis;
            ### Probably store it in a list. Ya know, for standard deviations and such
            ### Find the Binder cumulant
            # For A
            thisbinderzA = 1 - (mz4A/(3.*mz2A*mz2A))
            bindervalueszcompA.append(thisbinderzA)
            # For B
            thisbinderzB = 1 - (mz4B/(3.*mz2B*mz2B))
            bindervalueszcompB.append(thisbinderzB)
            
            ### Do a fitting of the results. To a line, a quadr. func., a pol. of the third degree. Check this. OFS suggested quadr
        if quadraticnotcubicfit==0:
            # For A
            bindervalueszcompA = array(bindervalueszcompA) 
            fitvectorzcompA = polyfit(betasA, bindervalueszcompA, 1) # Fits the function points to a quadratic polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1];
            # For B
            bindervalueszcompB = array(bindervalueszcompB) 
            fitvectorzcompB = polyfit(betasB, bindervalueszcompB, 1) # Fits the function points to a quadratic polynomial
            azB = fitvectorzcompB[0]; bzB = fitvectorzcompB[1]; 
                        
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = giveline(azA, bzA, fbetas) # Should this really be quadratic? I look at a small interval...
            fB = giveline(azB, bzB, fbetas)
        if quadraticnotcubicfit==1:
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
        if quadraticnotcubicfit==2:
            # For A
            bindervalueszcompA = array(bindervalueszcompA) 
            fitvectorzcompA = polyfit(betasA, bindervalueszcompA, 3) # Fits the function points to a quadratic polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]; dzA = fitvectorzcompA[3];
            # For B
            bindervalueszcompB = array(bindervalueszcompB) 
            fitvectorzcompB = polyfit(betasB, bindervalueszcompB, 3) # Fits the function points to a quadratic polynomial
            azB = fitvectorzcompB[0]; bzB = fitvectorzcompB[1]; czB = fitvectorzcompB[2]; dzB = fitvectorzcompB[3];
            
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = givecubic(azA, bzA, czA, dzA, fbetas) # Should this really be quadratic? I look at a small interval...
            fB = givecubic(azB, bzB, czB, dzB, fbetas)
        crossingbeta, crossingU = findintersection(fbetas, fA, fB)
        intersectionbeta_av += crossingbeta
        intersectionbetas[k] = crossingbeta
        intersectionU_av += crossingU
        intersectionUs[k] = crossingU
        
        if plotter==0:
            figure(figsize=(6,5))
            plot(fbetas, fA)
            hold('on')
            plot(fbetas, fB)
            plot(betasA, bindervalueszcompA, 'ro')
            plot(betasB, bindervalueszcompB, 'bo')
            title(r'Crossings between Binder cumulant graphs L and L+%i'%deltaL, fontsize=14)
            xlabel(r'$\beta$', fontsize=20)
            ylabel(r'$U_{L,z}$', fontsize=20)
            tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
            show()
            plotter = 1
        
    intersectionbeta_av = intersectionbeta_av/no_of_BootstraprunsA
    intersectionbeta_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        intersectionbeta_stddev += (intersectionbeta_av-intersectionbetas[i])*(intersectionbeta_av-intersectionbetas[i])
    intersectionbeta_stddev = sqrt(intersectionbeta_stddev/(no_of_BootstraprunsA*(no_of_BootstraprunsA-1)))
    
    intersectionU_av = intersectionU_av/no_of_BootstraprunsA
    intersectionU_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        intersectionU_stddev += (intersectionU_av-intersectionUs[i])*(intersectionU_av-intersectionUs[i])
    intersectionU_stddev = sqrt(intersectionU_stddev/(no_of_BootstraprunsA*(no_of_BootstraprunsA-1)))
    
    return intersectionbeta_av, intersectionbeta_stddev, intersectionU_av, intersectionU_stddev


#################### This is where we give our input #######################
LA = 6
LB = 8
LC = 10
LD = 12
LE = 14

'''
# Jensen's parameters  ## Please keep and comment out when using other param.s. Saving a lot of bother...
filenameA = "fcc6x6x6yopen_beta0p775to0p779_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameB = "fcc8x8x8yopen_beta0p775to0p779_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameC = "fcc10x10x10yopen_beta0p775to0p779_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameD = "fcc12x12x12yopen_beta0p775to0p779_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameE = "fcc14x14x14yopen_beta0p775to0p779_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
'''
# Li's parameters
filenameA = "fcc6x6x6yopen_beta0p7to0p82_Nbeta50_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameB = "fcc8x8x8yopen_beta0p7to0p82_Nbeta50_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameC = "fcc10x10x10yopen_beta0p7to0p82_Nbeta50_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameD = "fcc12x12x12yopen_beta0p7to0p82_Nbeta50_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameE = "fcc14x14x14yopen_beta0p7to0p82_Nbeta50_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"
filenameF = "fcc16x16x16yopen_beta0p7to0p82_Nbeta50_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_binavgs.txt"

##################                                   #######################
bool16 = 0 #0 - we have only up to L=14 / anythong else- we have L=16 too

quadraticnotcubicfit = 2 # 0 gives a line fit to the Binder cumulants, 2 gives a quadratic fit, 3 gives cubic
toplotornot = 0 # 0, don't plot, anything else leads to a plot

# If we only want to look at a small temperature interval
cutit = 1 # 0: Find a fit for the whole beta range; 1: Find a fit for a small interval
cutlower = 0.76
cutupper = 0.8

deltaL = 2 # Difference in the size of the systems for the graphs we use for the intersection

ib_avA, ib_stdvA, iU_avA, iU_stdvA = extract_crossing(filenameA, filenameB, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL)
ib_avB, ib_stdvB, iU_avB, iU_stdvB = extract_crossing(filenameB, filenameC, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL)
ib_avC, ib_stdvC, iU_avC, iU_stdvC = extract_crossing(filenameC, filenameD, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL)
ib_avD, ib_stdvD, iU_avD, iU_stdvD = extract_crossing(filenameD, filenameE, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL)
if bool16!=0:
    ib_avE, ib_stdvE, iU_avE, iU_stdvE = extract_crossing(filenameE, filenameF, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL)

crossingbetas      = [ib_avA, ib_avB, ib_avC, ib_avD]
crossingbetas_stdv = [ib_stdvA, ib_stdvB, ib_stdvC, ib_stdvD]
crossingUs         = [iU_avA, iU_avB, iU_avC, iU_avD]
crossingUs_stdv    = [iU_stdvA, iU_stdvB, iU_stdvC, iU_stdvD]
vec1dL             = [1./LA, 1./LB, 1./LC, 1./LD]

if bool16!=0:
    crossingbetas.append(ib_avE)
    crossingbetas_stdv.append(ib_stdvE)
    crossingUs.append(iU_avE)
    crossingUs_stdv.append(iU_stdvE)
    vec1dL.append(1./LE)

temps = []
for betaval in crossingbetas:
    temp = 11.6045221/betaval
    temps.append(temp)

crossingbetas      = array(crossingbetas)
crossingbetas_stdv = array(crossingbetas_stdv)
crossingUs         = array(crossingUs)
crossingUs_stdv    = array(crossingUs_stdv)
vec1dL             = array(vec1dL)
temps              = array(temps)

print "Crossing points of graphs"
for betapoint in crossingbetas:
    print betapoint

'''
print "Min crossingbetas", min(crossingbetas)
print "Max crossingbetas", max(crossingbetas)

print "Min crossingbetas*0.98", min(crossingbetas)*0.999
print "Max crossingbetas*1.02", max(crossingbetas)*1.001
'''

maxbeta_inplot = (max(vec1dL)*1.02)

# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL, crossingbetas, yerr=crossingbetas_stdv, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$\beta$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()

# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL, crossingUs, yerr=crossingUs_stdv, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, maxbeta_inplot, (min(crossingUs)-max(crossingUs_stdv))*0.999, (max(crossingUs)+max(crossingUs_stdv))*1.001])
show()


'''
# log-log plot
figure(figsize=(6,5))
#errorbar(vec1dL, crossingbetas, yerr=crossingbetas_stdv, fmt=None, capsize=2)
loglog(vec1dL, crossingbetas, basex=2)
#semilogy(vec1dL, crossingbetas)
#semilogx(vec1dL, crossingbetas)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$\beta$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()

# log-log plot
figure(figsize=(6,5))
#errorbar(vec1dL, crossingbetas, yerr=crossingbetas_stdv, fmt=None, capsize=2)
#loglog(vec1dL, crossingbetas, basex=2)
semilogy(vec1dL, crossingbetas)
#semilogx(vec1dL, crossingbetas)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$\beta$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()
'''
'''
# log-log plot
figure(figsize=(6,5))
#errorbar(vec1dL, crossingbetas, yerr=crossingbetas_stdv, fmt=None, capsize=2)
#loglog(vec1dL, crossingbetas, basex=2)
#semilogy(vec1dL, crossingbetas)
semilogx(vec1dL, crossingbetas)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$\beta$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()
'''


# Finding the exponent a la the log-log plot
# From www.physics.pomona.edu/sixideas/old/labs/LRM/LR05.pdf
#n = (log(crossingbetas[len(crossingbetas)-1])-log(crossingbetas[0]))/(log(vec1dL[len(vec1dL)-1])-log(vec1dL[0]))
n = (log(temps[0])-log(temps[len(crossingbetas)-1]))/(log(vec1dL[0])-log(vec1dL[len(vec1dL)-1]))
endpointbeta = log(temps[len(crossingbetas)-1])
endpointinverselength = log(vec1dL[len(vec1dL)-1])

print " "
print "slope of crossing point graph, log-log =", n
#print "log(10)=", log(10), "log(2)=", log(2), "log(e)=", log(2.7182818284590452353602874713527)
k = exp(endpointbeta-(n*endpointinverselength))
print "k =", k
print " "


# Trying to find where it crosses the y-axis
fitLs = linspace(0,maxbeta_inplot,1000)
# Different degrees
# First degree
fq0 = polyfit(vec1dL, crossingbetas, 1)
a1 = fq0[0]; b1 = fq0[1];
fit_line      = giveline(a1,b1,fitLs)
# Second degree
fq1 = polyfit(vec1dL, crossingbetas, 2)
a2 = fq1[0]; b2 = fq1[1]; c2 = fq1[2]
fit_quadratic = givequadratic(a2,b2,c2,fitLs)
# Third degree
fq2 = polyfit(vec1dL, crossingbetas, 3)
a3 = fq2[0]; b3 = fq2[1]; c3 = fq2[2]; d3 = fq2[3]
fit_cubic     = givecubic(a3,b3,c3,d3,fitLs)
'''
# Fourth degree
fq3 = polyfit(vec1dL, crossingbetas, 4)
a4 = fq3[0]; b4 = fq3[1]; c4 = fq3[2]; d4 = fq3[3]; e4 = fq3[4]
fit_quadruple = givequadruple(a4,b4,c4,d4,e4,fitLs)
'''
#Power law
#fit_powerlaw = givepowerlaw(k,n,fitLs)

print "Estimated crossings (beta):"
print "Line fit:", b1
print "Quadratic polynomial:", c2
print "Cubic polynomial:", d3
#print "Quadratic polynomial:", e4
#print "Power law:", fit_powerlaw[0]

figure(figsize=(6,5))
errorbar(vec1dL, crossingbetas, yerr=crossingbetas_stdv, fmt=None, capsize=2)
hold('on')
plot(fitLs, fit_line, label='Line')
plot(fitLs, fit_quadratic, label='Quadratic')
plot(fitLs, fit_cubic, label='Cubic')
#plot(fitLs, fit_quadruple, label='Quadruple')
#plot(fitLs, fit_powerlaw, label='Power law')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$\beta$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower center")
show()

# Do something like crossingbetas+errorbars, crossingbetas-errorbars etc to get errorbars

#################### BOOTSTRAP L, L+4######################################

quadraticnotcubicfit = 2 # 0 gives a line fit to the Binder cumulants, 2 gives a quadratic fit, 3 gives cubic
toplotornot = 1 # 0, don't plot, anything else leads to a plot
# A-6, B-8, C-10, D-12, E-14, F-16
# 6-10-14, A-C-E
# 8-12-16, B-D-F


deltaL = 4 # Difference in the size of the systems for the graphs we use for the intersection

ib_avA4, ib_stdvA4, iU_avA4, iU_stdvA4 = extract_crossing(filenameA, filenameC, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL)
ib_avB4, ib_stdvB4, iU_avB4, iU_stdvB4 = extract_crossing(filenameB, filenameD, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL)
ib_avC4, ib_stdvC4, iU_avC4, iU_stdvC4 = extract_crossing(filenameC, filenameE, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL)

if bool16!=0:
    ib_avD4, ib_stdvD4, iU_avD4, iU_stdvD4 = extract_crossing(filenameD, filenameF, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL)

crossingbetas4      = [ib_avA4, ib_avB4, ib_avC4]
crossingbetas_stdv4 = [ib_stdvA4, ib_stdvB4, ib_stdvC4]
crossingUs4         = [iU_avA4, iU_avB4, iU_avC4]
crossingUs_stdv4    = [iU_stdvA4, iU_stdvB4, iU_stdvC4]
vec1dL4             = [1./LA, 1./LB, 1./LC]

if bool16!=0:
    crossingbetas4.append(ib_avD4)
    crossingbetas_stdv4.append(ib_stdvD4)
    crossingUs4.append(iU_avD4)
    crossingUs_stdv4.append(iU_stdvD4)
    vec1dL4.append(1./LD)

temps4 = []
for betaval in crossingbetas4:
    temp = 11.6045221/betaval
    temps4.append(temp)

crossingbetas4      = array(crossingbetas4)
crossingbetas_stdv4 = array(crossingbetas_stdv4)
crossingUs4         = array(crossingUs4)
crossingUs_stdv4    = array(crossingUs_stdv4)
vec1dL4             = array(vec1dL4)
temps4              = array(temps4)

print "Crossing points of graphs"
for betapoint in crossingbetas4:
    print betapoint

'''
print "Min crossingbetas", min(crossingbetas)
print "Max crossingbetas", max(crossingbetas)

print "Min crossingbetas*0.98", min(crossingbetas)*0.999
print "Max crossingbetas*1.02", max(crossingbetas)*1.001
'''

maxbeta_inplot = (max(vec1dL4)*1.02)

# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL4, crossingbetas4, yerr=crossingbetas_stdv4, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+4 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$\beta$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, maxbeta_inplot, min(crossingbetas4)*0.999, max(crossingbetas4)*1.001])
show()

# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL4, crossingUs4, yerr=crossingUs_stdv4, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+4 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, maxbeta_inplot, (min(crossingUs4)-max(crossingUs_stdv4))*0.999, (max(crossingUs4)+max(crossingUs_stdv4))*1.001])
show()

# Trying to find where it crosses the y-axis
fitLs = linspace(0,maxbeta_inplot,1000)
# Different degrees
# First degree
fq0 = polyfit(vec1dL4, crossingbetas4, 1)
a1 = fq0[0]; b1 = fq0[1];
fit_line      = giveline(a1,b1,fitLs)
# Second degree
fq1 = polyfit(vec1dL4, crossingbetas4, 2)
a2 = fq1[0]; b2 = fq1[1]; c2 = fq1[2]
fit_quadratic = givequadratic(a2,b2,c2,fitLs)
# Third degree
fq2 = polyfit(vec1dL4, crossingbetas4, 3)
a3 = fq2[0]; b3 = fq2[1]; c3 = fq2[2]; d3 = fq2[3]
fit_cubic     = givecubic(a3,b3,c3,d3,fitLs)
'''
# Fourth degree
fq3 = polyfit(vec1dL, crossingbetas, 4)
a4 = fq3[0]; b4 = fq3[1]; c4 = fq3[2]; d4 = fq3[3]; e4 = fq3[4]
fit_quadruple = givequadruple(a4,b4,c4,d4,e4,fitLs)
'''
#Power law
#fit_powerlaw = givepowerlaw(k,n,fitLs)

print "Estimated crossings (beta), from intersections of graphs with L and L+4:"
print "Line fit:", b1
print "Quadratic polynomial:", c2
print "Cubic polynomial:", d3
#print "Quadratic polynomial:", e4
#print "Power law:", fit_powerlaw[0]

figure(figsize=(6,5))
errorbar(vec1dL4, crossingbetas4, yerr=crossingbetas_stdv4, fmt=None, capsize=2)
hold('on')
plot(fitLs, fit_line, label='Line')
plot(fitLs, fit_quadratic, label='Quadratic')
plot(fitLs, fit_cubic, label='Cubic')
#plot(fitLs, fit_quadruple, label='Quadruple')
#plot(fitLs, fit_powerlaw, label='Power law')
title(r'Crossings between Binder cumulant graphs L and L+4 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$\beta$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower center")
show()



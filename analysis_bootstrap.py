from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys
from random import randint

# Our system sizes # I guess I have too many
LA = 6
LB = 8
LC = 10
LD = 12
LE = 12
LF = 14
#LG = 20

# ----------------------------------------------------A-------------------------------------------------#
# Open the file for reading # Give the correct address

def get_Bootstrapresults(filename):
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
    
    # Should I work with lists or arrays?# Make an array for the results from the bootstrap. Should I do it here, or?
    
    avectorzA = []; bvectorzA = []; cvectorzA = []
    azA_av = 0; bzA_av = 0; czA_av = 0;
    
    ### Oy, vey, all the loops... # Need one for several bootstrap runs too.
    ### for different numbers of bins for different beta
    for k in range(0, no_of_Bootstrapruns): # Want to repeat the random selection of bins a number of times
        bindervalueszcomp = []
        for i in range(0, no_of_betasA): # For doing it for every temperature
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
            ### Find the Binder cumulant
            thisbinderz = 1 - (mz4/(3.*mz2*mz2))
            bindervalueszcomp.append(thisbinderz) 
            ### Find other stuff as well 
            ### Do a fitting of the results. To a line, a quadr. func., a pol. of the third degree. Check this. OFS suggested quadr
        bindervalueszcomp = array(bindervalueszcomp) 
        fitvectorzcomp = polyfit(betasA, bindervalueszcomp, 2) # Fits the function points to a quadratic polynomial
        ### Store the coeffs of the fitting
        az = fitvectorzcomp[0]; bz = fitvectorzcomp[1]; cz = fitvectorzcomp[2]
        azA_av += az; bzA_av += bz; czA_av += cz;
        avectorzA.append(az)
        bvectorzA.append(bz)
        cvectorzA.append(cz)
    azA_av = azA_av/no_of_Bootstrapruns; bzA_av = bzA_av/no_of_Bootstrapruns; czA_av = czA_av/no_of_Bootstrapruns;
    
    azA_stdv = 0; bzA_stdv = 0; czA_stdv = 0;
    for i in range(0, no_of_Bootstrapruns):   # Is this the correct method?
        azA_stdv += (avectorzA[i]-azA_av)*(avectorzA[i]-azA_av)
        bzA_stdv += (bvectorzA[i]-bzA_av)*(bvectorzA[i]-bzA_av)
        czA_stdv += (cvectorzA[i]-czA_av)*(cvectorzA[i]-czA_av)
    azA_stdv = sqrt(azA_stdv/(no_of_Bootstrapruns*(no_of_Bootstrapruns-1)))
    bzA_stdv = sqrt(bzA_stdv/(no_of_Bootstrapruns*(no_of_Bootstrapruns-1)))
    czA_stdv = sqrt(czA_stdv/(no_of_Bootstrapruns*(no_of_Bootstrapruns-1)))
    betamin = betasA[0]
    betamax = betasA[len(betasA)-1]
    return betamin, betamax, azA_av, bzA_av, czA_av, azA_av, bzA_av, czA_av # Add more output as implement more...
    
def givequadratic(a, b, c, x):
    return a*x**2+b*x+c

def findintersection(beta, f1, f2):
    diff = f1 - f2
    for i in range(len(diff)):
        if (diff[i] == 0) or (diff[i]*diff[i+1] < 0):
            return beta[i]

filenameA = "test_binavgs.txt"
filenameB = "test_binavgs.txt"
filenameC = "test_binavgs.txt"
filenameD = "test_binavgs.txt"
betaminA, betamaxA, azA_av, bzA_av, czA_av, azA_av, bzA_av, czA_av = get_Bootstrapresults(filenameA)
betaminB, betamaxB, azB_av, bzB_av, czB_av, azB_av, bzB_av, czB_av = get_Bootstrapresults(filenameA)
betaminC, betamaxC, azC_av, bzC_av, czC_av, azC_av, bzC_av, czC_av = get_Bootstrapresults(filenameA)
betaminD, betamaxD, azD_av, bzD_av, czD_av, azD_av, bzD_av, czD_av = get_Bootstrapresults(filenameA)

res = 1000 # Resolution of our arrays
betasA = linspace(betaminA, betamaxA, res)
betasB = linspace(betaminB, betamaxB, res)
betasC = linspace(betaminC, betamaxC, res)
betasD = linspace(betaminD, betamaxD, res)
fA = givequadratic(azA_av, bzA_av, czA_av, betasA)
fB = givequadratic(azB_av, bzB_av, czB_av, betasB)
fC = givequadratic(azC_av, bzC_av, czC_av, betasC)
fD = givequadratic(azD_av, bzD_av, czD_av, betasD)

# But how do I deal with the uncertainties in a,b and c? That must enter in the final answer, somehow
testf = zeros(res)+0.4
testbeta = findintersection(betasA, fA, testf)
print testbeta

# For testing
figure()
plot(betasA, fA)
hold('on')
plot(betasA, testf)
show()

'''
# For testing
figure()
plot(betasA, fA)
hold('on')
plot(betasB, fB)
plot(betasC, fC)
plot(betasD, fD)
show()
'''

'''
for i in range(0, len(binofbeta_lastindexA)):
    print "binofbeta_lastindexA[", i, "] = " , binofbeta_lastindexA[i]
'''

# Should I have all of the above as a kind of function? Copy-pasting is a bit impractical. But how to do it with a file...
# I wouldn't want to overload myself with data...
# Could I take infile as input?

# Loop over beta first?
# Probably have a few functions to make calls to.

# Need to do this for several system sizes

# THEN find the intersections of the graphs


# Rough draft on plotting
'''
figure(figsize=(6,5))
errorbar(betasA, ULAy, yerr=ULAdelta_y, capsize=2, label='L=%i'%LA)
hold('on')
errorbar(betasB, ULBy, yerr=ULBdelta_y, capsize=2, label='L=%i'%LB)
errorbar(betasC, ULCy, yerr=ULCdelta_y, capsize=2, label='L=%i'%LC)
errorbar(betasD, ULDy, yerr=ULDdelta_y, capsize=2, label='L=%i'%LD)
#errorbar(betasE, ULEy, yerr=ULEdelta_y, capsize=2, label='L=%i'%LE)
#errorbar(betasF, ULFy, yerr=ULFdelta_y, capsize=2, label='L=%i'%LF)
#errorbar(betasG, ULGy, yerr=ULGdelta_y, capsize=2, label='L=%i'%LG)
title(r'The Binder cumulant vs $\beta$ for different L', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$U_{L,y}$', fontsize=20)
legend(loc="lower left")
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
show()
'''

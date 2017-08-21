from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Our system sizes # I guess I have too many
LA = 6
LB = 8
LC = 16
LD = 20
#LE = 14
#LF = 16
#LG = 20

# ----------------------------------------------------A-------------------------------------------------#
# Open the file for reading # Give the correct address
infile = open("fcc6x6x6yopen_beta0p1to5_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000_bins100_seed79_latticeseed21_slowcool_mxyzq2pi010.txt", "r")

# Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
betasA            = []          # List of beta values
mxsq_avsA         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmxsqsA         = []          # Average m over all bins and mcsteps, listed according to beta value
mxquadravsA       = []
mxquadravs_stdA   = []
mysq_avsA         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmysqsA         = []          # Average m over all bins and mcsteps, listed according to beta value
myquadravsA       = []
myquadravs_stdA   = []
mzsq_avsA         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmzsqsA         = []          # Average m over all bins and mcsteps, listed according to beta value
mzquadravsA       = []
mzquadravs_stdA   = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        beta = float(words[0])
        betasA.append(beta)
        # mxsq_avs                    # <m_x^2>
        mxsq_av = float(words[1])
        mxsq_avsA.append(mxsq_av)
        # stdmxsqsA
        stdm = float(words[2])
        stdmxsqsA.append(stdm)
        # mysq_avs                    # <m_y^2>
        mysq_av = float(words[3])
        mysq_avsA.append(mysq_av)
        # stdmxsqsA
        stdm = float(words[4])
        stdmysqsA.append(stdm)
        # mzsq_avs                    # <m_z^2>
        mzsq_av = float(words[5])
        mzsq_avsA.append(mzsq_av)
        # stdmzsqsA
        stdm = float(words[6])
        stdmzsqsA.append(stdm)
        # mxquadravs20                 # <m_x^4>
        mq = float(words[7])
        mxquadravsA.append(mq)
        # mxquadravs_std20
        mq_std = float(words[8])
        mxquadravs_stdA.append(mq_std)
        # myquadravs20                 # <m_y^4>
        mq = float(words[9])
        myquadravsA.append(mq)
        # mxquadravs_stdA
        mq_std = float(words[10])
        myquadravs_stdA.append(mq_std)
        # mzquadravs20                 # <m_z^4>
        mq = float(words[11])
        mzquadravsA.append(mq)
        # mzquadravs_std20
        mq_std = float(words[12])
        mzquadravs_stdA.append(mq_std)

# Remember to close the file
infile.close()


# ---------------------------------------------------------B--------------------------------------------------
# Open the file for reading # Give the correct address
infile = open("fcc8x8x8yopen_beta0p1to5_Nbeta50_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000_bins100_seed79_latticeseed21_slowcool_mxyzq2pi010.txt", "r")

# Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
betasB            = []          # List of beta values
mxsq_avsB         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmxsqsB         = []          # Average m over all bins and mcsteps, listed according to beta value
mxquadravsB       = []
mxquadravs_stdB   = []
mysq_avsB         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmysqsB         = []          # Average m over all bins and mcsteps, listed according to beta value
myquadravsB       = []
myquadravs_stdB   = []
mzsq_avsB         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmzsqsB         = []          # Average m over all bins and mcsteps, listed according to beta value
mzquadravsB       = []
mzquadravs_stdB   = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        beta = float(words[0])
        betasB.append(beta)
        # mxsq_avs                    # <m_x^2>
        mxsq_av = float(words[1])
        mxsq_avsB.append(mxsq_av)
        # stdmxsqsB
        stdm = float(words[2])
        stdmxsqsB.append(stdm)
        # mysq_avs                    # <m_y^2>
        mysq_av = float(words[3])
        mysq_avsB.append(mysq_av)
        # stdmxsqsB
        stdm = float(words[4])
        stdmysqsB.append(stdm)
        # mzsq_avs                    # <m_z^2>
        mzsq_av = float(words[5])
        mzsq_avsB.append(mzsq_av)
        # stdmzsqsB
        stdm = float(words[6])
        stdmzsqsB.append(stdm)
        # mxquadravs20                 # <m_x^4>
        mq = float(words[7])
        mxquadravsB.append(mq)
        # mxquadravs_std20
        mq_std = float(words[8])
        mxquadravs_stdB.append(mq_std)
        # myquadravs20                 # <m_y^4>
        mq = float(words[9])
        myquadravsB.append(mq)
        # mxquadravs_stdA
        mq_std = float(words[10])
        myquadravs_stdB.append(mq_std)
        # mzquadravs20                 # <m_z^4>
        mq = float(words[11])
        mzquadravsB.append(mq)
        # mzquadravs_std20
        mq_std = float(words[12])
        mzquadravs_stdB.append(mq_std)

# Remember to close the file
infile.close()

'''
# --------------------------------------------------------C----------------------------------------------
# Open the file for reading # Give the correct address
infile = open("test_mxyzq2pi010.txt", "r")

# Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
betasC            = []          # List of beta values
mxsq_avsC         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmxsqsC         = []          # Average m over all bins and mcsteps, listed according to beta value
mxquadravsC       = []
mxquadravs_stdC   = []
mysq_avsC         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmysqsC         = []          # Average m over all bins and mcsteps, listed according to beta value
myquadravsC       = []
myquadravs_stdC   = []
mzsq_avsC         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmzsqsC         = []          # Average m over all bins and mcsteps, listed according to beta value
mzquadravsC       = []
mzquadravs_stdC   = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        beta = float(words[0])
        betasC.append(beta)
        # mxsq_avs                    # <m_x^2>
        mxsq_av = float(words[1])
        mxsq_avsC.append(mxsq_av)
        # stdmxsqsA
        stdm = float(words[2])
        stdmxsqsC.append(stdm)
        # mysq_avs                    # <m_y^2>
        mysq_av = float(words[3])
        mysq_avsC.append(mysq_av)
        # stdmxsqsA
        stdm = float(words[4])
        stdmysqsC.append(stdm)
        # mzsq_avs                    # <m_z^2>
        mzsq_av = float(words[5])
        mzsq_avsC.append(mzsq_av)
        # stdmzsqsA
        stdm = float(words[6])
        stdmzsqsC.append(stdm)
        # mxquadravs20                 # <m_x^4>
        mq = float(words[7])
        mxquadravsC.append(mq)
        # mxquadravs_std20
        mq_std = float(words[8])
        mxquadravs_stdC.append(mq_std)
        # myquadravs20                 # <m_y^4>
        mq = float(words[9])
        myquadravsC.append(mq)
        # mxquadravs_stdA
        mq_std = float(words[10])
        myquadravs_stdC.append(mq_std)
        # mzquadravs20                 # <m_z^4>
        mq = float(words[11])
        mzquadravsC.append(mq)
        # mzquadravs_std20
        mq_std = float(words[12])
        mzquadravs_stdC.append(mq_std)

# Remember to close the file
infile.close()

# --------------------------------------------D---------------------------------
# Open the file for reading # Give the correct address
infile = open("test_mxyzq2pi010.txt", "r")

# Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
betasD            = []          # List of beta values
mxsq_avsD         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmxsqsD         = []          # Average m over all bins and mcsteps, listed according to beta value
mxquadravsD       = []
mxquadravs_stdD   = []
mysq_avsD         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmysqsD         = []          # Average m over all bins and mcsteps, listed according to beta value
myquadravsD       = []
myquadravs_stdD   = []
mzsq_avsD         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmzsqsD         = []          # Average m over all bins and mcsteps, listed according to beta value
mzquadravsD       = []
mzquadravs_stdD   = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        beta = float(words[0])
        betasD.append(beta)
        # mxsq_avs                    # <m_x^2>
        mxsq_av = float(words[1])
        mxsq_avsD.append(mxsq_av)
        # stdmxsqsA
        stdm = float(words[2])
        stdmxsqsD.append(stdm)
        # mysq_avs                    # <m_y^2>
        mysq_av = float(words[3])
        mysq_avsD.append(mysq_av)
        # stdmxsqsA
        stdm = float(words[4])
        stdmysqsD.append(stdm)
        # mzsq_avs                    # <m_z^2>
        mzsq_av = float(words[5])
        mzsq_avsD.append(mzsq_av)
        # stdmzsqsA
        stdm = float(words[6])
        stdmzsqsD.append(stdm)
        # mxquadravs20                 # <m_x^4>
        mq = float(words[7])
        mxquadravsD.append(mq)
        # mxquadravs_std20
        mq_std = float(words[8])
        mxquadravs_stdD.append(mq_std)
        # myquadravs20                 # <m_y^4>
        mq = float(words[9])
        myquadravsD.append(mq)
        # mxquadravs_stdA
        mq_std = float(words[10])
        myquadravs_stdD.append(mq_std)
        # mzquadravs20                 # <m_z^4>
        mq = float(words[11])
        mzquadravsD.append(mq)
        # mzquadravs_std20
        mq_std = float(words[12])
        mzquadravs_stdD.append(mq_std)

# Remember to close the file
infile.close()

# --------------------------------------------E------------------------------
# Open the file for reading # Give the correct address
infile = open("test_mxyzq2pi010.txt", "r")

# Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
betasE            = []          # List of beta values
mxsq_avsE         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmxsqsE         = []          # Average m over all bins and mcsteps, listed according to beta value
mxquadravsE       = []
mxquadravs_stdE   = []
mysq_avsE         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmysqsE         = []          # Average m over all bins and mcsteps, listed according to beta value
myquadravsE       = []
myquadravs_stdE   = []
mzsq_avsE         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmzsqsE         = []          # Average m over all bins and mcsteps, listed according to beta value
mzquadravsE       = []
mzquadravs_stdE   = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        beta = float(words[0])
        betasE.append(beta)
        # mxsq_avs                    # <m_x^2>
        mxsq_av = float(words[1])
        mxsq_avsE.append(mxsq_av)
        # stdmxsqsA
        stdm = float(words[2])
        stdmxsqsE.append(stdm)
        # mysq_avs                    # <m_y^2>
        mysq_av = float(words[3])
        mysq_avsE.append(mysq_av)
        # stdmxsqsA
        stdm = float(words[4])
        stdmysqsE.append(stdm)
        # mzsq_avs                    # <m_z^2>
        mzsq_av = float(words[5])
        mzsq_avsE.append(mzsq_av)
        # stdmzsqsA
        stdm = float(words[6])
        stdmzsqsE.append(stdm)
        # mxquadravs20                 # <m_x^4>
        mq = float(words[7])
        mxquadravsE.append(mq)
        # mxquadravs_std20
        mq_std = float(words[8])
        mxquadravs_stdE.append(mq_std)
        # myquadravs20                 # <m_y^4>
        mq = float(words[9])
        myquadravsE.append(mq)
        # mxquadravs_stdA
        mq_std = float(words[10])
        myquadravs_stdE.append(mq_std)
        # mzquadravs20                 # <m_z^4>
        mq = float(words[11])
        mzquadravsE.append(mq)
        # mzquadravs_std20
        mq_std = float(words[12])
        mzquadravs_stdE.append(mq_std)

# Remember to close the file
infile.close()



# --------------------------------------------G---------------------------------------------
# Open the file for reading # Give the correct address
infile = open("test_mxyzq2pi010.txt", "r")

# Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
betasG            = []          # List of beta values
mxsq_avsG         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmxsqsG         = []          # Average m over all bins and mcsteps, listed according to beta value
mxquadravsG       = []
mxquadravs_stdG   = []
mysq_avsG         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmysqsG         = []          # Average m over all bins and mcsteps, listed according to beta value
myquadravsG       = []
myquadravs_stdG   = []
mzsq_avsG         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmzsqsG         = []          # Average m over all bins and mcsteps, listed according to beta value
mzquadravsG       = []
mzquadravs_stdG   = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        beta = float(words[0])
        betasG.append(beta)
        # mxsq_avs                    # <m_x^2>
        mxsq_av = float(words[1])
        mxsq_avsG.append(mxsq_av)
        # stdmxsqsA
        stdm = float(words[2])
        stdmxsqsG.append(stdm)
        # mysq_avs                    # <m_y^2>
        mysq_av = float(words[3])
        mysq_avsG.append(mysq_av)
        # stdmxsqsA
        stdm = float(words[4])
        stdmysqsG.append(stdm)
        # mzsq_avs                    # <m_z^2>
        mzsq_av = float(words[5])
        mzsq_avsG.append(mzsq_av)
        # stdmzsqsA
        stdm = float(words[6])
        stdmzsqsG.append(stdm)
        # mxquadravs20                 # <m_x^4>
        mq = float(words[7])
        mxquadravsG.append(mq)
        # mxquadravs_std20
        mq_std = float(words[8])
        mxquadravs_stdG.append(mq_std)
        # myquadravs20                 # <m_y^4>
        mq = float(words[9])
        myquadravsG.append(mq)
        # mxquadravs_stdA
        mq_std = float(words[10])
        myquadravs_stdG.append(mq_std)
        # mzquadravs20                 # <m_z^4>
        mq = float(words[11])
        mzquadravsG.append(mq)
        # mzquadravs_std20
        mq_std = float(words[12])
        mzquadravs_stdG.append(mq_std)

# Remember to close the file
infile.close()

# ------------------------------------------------F-------------------------------------
# Open the file for reading # Give the correct address
infile = open("test_mxyzq2pi010.txt", "r")

# Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
betasF            = []          # List of beta values
mxsq_avsF         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmxsqsF         = []          # Average m over all bins and mcsteps, listed according to beta value
mxquadravsF       = []
mxquadravs_stdF   = []
mysq_avsF         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmysqsF         = []          # Average m over all bins and mcsteps, listed according to beta value
myquadravsF       = []
myquadravs_stdF   = []
mzsq_avsF         = []          # Average m over all bins and mcsteps, listed according to beta value
stdmzsqsF         = []          # Average m over all bins and mcsteps, listed according to beta value
mzquadravsF       = []
mzquadravs_stdF   = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        beta = float(words[0])
        betasF.append(beta)
        # mxsq_avs                    # <m_x^2>
        mxsq_av = float(words[1])
        mxsq_avsF.append(mxsq_av)
        # stdmxsqsA
        stdm = float(words[2])
        stdmxsqsF.append(stdm)
        # mysq_avs                    # <m_y^2>
        mysq_av = float(words[3])
        mysq_avsF.append(mysq_av)
        # stdmxsqsA
        stdm = float(words[4])
        stdmysqsF.append(stdm)
        # mzsq_avs                    # <m_z^2>
        mzsq_av = float(words[5])
        mzsq_avsF.append(mzsq_av)
        # stdmzsqsA
        stdm = float(words[6])
        stdmzsqsF.append(stdm)
        # mxquadravs20                 # <m_x^4>
        mq = float(words[7])
        mxquadravsF.append(mq)
        # mxquadravs_std20
        mq_std = float(words[8])
        mxquadravs_stdF.append(mq_std)
        # myquadravs20                 # <m_y^4>
        mq = float(words[9])
        myquadravsF.append(mq)
        # mxquadravs_stdA
        mq_std = float(words[10])
        myquadravs_stdF.append(mq_std)
        # mzquadravs20                 # <m_z^4>
        mq = float(words[11])
        mzquadravsF.append(mq)
        # mzquadravs_std20
        mq_std = float(words[12])
        mzquadravs_stdF.append(mq_std)

# Remember to close the file
infile.close()
'''

betasA            = array(betasA)
mxsq_avsA         = array(mxsq_avsA)
stdmxsqsA         = array(stdmxsqsA)
mxquadravsA       = array(mxquadravsA)
mxquadravs_stdA   = array(mxquadravs_stdA)
mysq_avsA         = array(mysq_avsA)
stdmysqsA         = array(stdmysqsA)
myquadravsA       = array(myquadravsA)
myquadravs_stdA   = array(myquadravs_stdA)
mzsq_avsA         = array(mzsq_avsA)
stdmzsqsA         = array(stdmzsqsA)
mzquadravsA       = array(mzquadravsA)
mzquadravs_stdA   = array(mzquadravs_stdA)


betasB            = array(betasB)
mxsq_avsB         = array(mxsq_avsB)
stdmxsqsB         = array(stdmxsqsB)
mxquadravsB       = array(mxquadravsB)
mxquadravs_stdB   = array(mxquadravs_stdB)
mysq_avsB         = array(mysq_avsB)
stdmysqsB         = array(stdmysqsB)
myquadravsB       = array(myquadravsB)
myquadravs_stdB   = array(myquadravs_stdB)
mzsq_avsB         = array(mzsq_avsB)
stdmzsqsB         = array(stdmzsqsB)
mzquadravsB       = array(mzquadravsB)
mzquadravs_stdB   = array(mzquadravs_stdB)

'''
betasC            = array(betasC)
mxsq_avsC         = array(mxsq_avsC)
stdmxsqsC         = array(stdmxsqsC)
mxquadravsC       = array(mxquadravsC)
mxquadravs_stdC   = array(mxquadravs_stdC)
mysq_avsC         = array(mysq_avsC)
stdmysqsC         = array(stdmysqsC)
myquadravsC       = array(myquadravsC)
myquadravs_stdC   = array(myquadravs_stdC)
mzsq_avsC         = array(mzsq_avsC)
stdmzsqsC         = array(stdmzsqsC)
mzquadravsC       = array(mzquadravsC)
mzquadravs_stdC   = array(mzquadravs_stdC)

betasD            = array(betasD)
mxsq_avsD         = array(mxsq_avsD)
stdmxsqsD         = array(stdmxsqsD)
mxquadravsD       = array(mxquadravsD)
mxquadravs_stdD   = array(mxquadravs_stdD)
mysq_avsD         = array(mysq_avsD)
stdmysqsD         = array(stdmysqsD)
myquadravsD       = array(myquadravsD)
myquadravs_stdD   = array(myquadravs_stdD)
mzsq_avsD         = array(mzsq_avsD)
stdmzsqsD         = array(stdmzsqsD)
mzquadravsD       = array(mzquadravsD)
mzquadravs_stdD   = array(mzquadravs_stdD)

betasE            = array(betasE)
mxsq_avsE         = array(mxsq_avsE)
stdmxsqsE         = array(stdmxsqsE)
mxquadravsE       = array(mxquadravsE)
mxquadravs_stdE   = array(mxquadravs_stdE)
mysq_avsE         = array(mysq_avsE)
stdmysqsE         = array(stdmysqsE)
myquadravsE       = array(myquadravsE)
myquadravs_stdE   = array(myquadravs_stdE)
mzsq_avsE         = array(mzsq_avsE)
stdmzsqsE         = array(stdmzsqsE)
mzquadravsE       = array(mzquadravsE)
mzquadravs_stdE   = array(mzquadravs_stdE)

betasF            = array(betasF)
mxsq_avsF         = array(mxsq_avsF)
stdmxsqsF         = array(stdmxsqsF)
mxquadravsF       = array(mxquadravsF)
mxquadravs_stdF   = array(mxquadravs_stdF)
mysq_avsF         = array(mysq_avsF)
stdmysqsF         = array(stdmysqsF)
myquadravsF       = array(myquadravsF)
myquadravs_stdF   = array(myquadravs_stdF)
mzsq_avsF         = array(mzsq_avsF)
stdmzsqsF         = array(stdmzsqsF)
mzquadravsF       = array(mzquadravsF)
mzquadravs_stdF   = array(mzquadravs_stdF)

betasG            = array(betasG)
mxsq_avsG         = array(mxsq_avsG)
stdmxsqsG         = array(stdmxsqsG)
mxquadravsG       = array(mxquadravsG)
mxquadravs_stdG   = array(mxquadravs_stdG)
mysq_avsG         = array(mysq_avsG)
stdmysqsG         = array(stdmysqsG)
myquadravsG       = array(myquadravsG)
myquadravs_stdG   = array(myquadravs_stdG)
mzsq_avsG         = array(mzsq_avsG)
stdmzsqsG         = array(stdmzsqsG)
mzquadravsG       = array(mzquadravsG)
mzquadravs_stdG   = array(mzquadravs_stdG)
'''

# I need to know how I should treat 3D-spin Binder cumulants

lengthA = len(betasA)
lengthB = len(betasB)
'''
lengthC = len(betasC)
lengthD = len(betasD)
lengthE = len(betasE)
lengthF = len(betasF)
lengthG = len(betasG)
'''

ULAx = zeros(lengthA)
ULAy = zeros(lengthA)
ULAz = zeros(lengthA)
ULBx = zeros(lengthB)
ULBy = zeros(lengthB)
ULBz = zeros(lengthB)
'''
ULCx = zeros(lengthC)
ULCy = zeros(lengthC)
ULCz = zeros(lengthC)
ULDx = zeros(lengthD)
ULDy = zeros(lengthD)
ULDz = zeros(lengthD)
ULEx = zeros(lengthE)
ULEy = zeros(lengthE)
ULEz = zeros(lengthE)
ULFx = zeros(lengthF)
ULFy = zeros(lengthF)
ULFz = zeros(lengthF)
ULGx = zeros(lengthG)
ULGy = zeros(lengthG)
ULGz = zeros(lengthG)
'''

ULAdelta_x = zeros(lengthA)
ULAdelta_y = zeros(lengthA)
ULAdelta_z = zeros(lengthA)
ULBdelta_x = zeros(lengthB)
ULBdelta_y = zeros(lengthB)
ULBdelta_z = zeros(lengthB)
'''
ULCdelta_x = zeros(lengthC)
ULCdelta_y = zeros(lengthC)
ULCdelta_z = zeros(lengthC)
ULDdelta_x = zeros(lengthD)
ULDdelta_y = zeros(lengthD)
ULDdelta_z = zeros(lengthD)
ULEdelta_x = zeros(lengthE)
ULEdelta_y = zeros(lengthE)
ULEdelta_z = zeros(lengthE)
ULFdelta_x = zeros(lengthF)
ULFdelta_y = zeros(lengthF)
ULFdelta_z = zeros(lengthF)
ULGdelta_x = zeros(lengthG)
ULGdelta_y = zeros(lengthG)
ULGdelta_z = zeros(lengthG)
'''

   
for i in range(lengthA):
    ULAx[i] = 1 - (mxquadravsA[i]/(3*mxsq_avsA[i]**2))
    ULAy[i] = 1 - (myquadravsA[i]/(3*mysq_avsA[i]**2))
    ULAz[i] = 1 - (mzquadravsA[i]/(3*mzsq_avsA[i]**2))
    #
    ULAdelta_x[i] = ULAx[i]*sqrt((mxquadravs_stdA[i]/mxquadravsA[i])**2+(stdmxsqsA[i]/mxsq_avsA[i])**2)
    ULAdelta_y[i] = ULAy[i]*sqrt((myquadravs_stdA[i]/myquadravsA[i])**2+(stdmysqsA[i]/mysq_avsA[i])**2)
    ULAdelta_z[i] = ULAz[i]*sqrt((mzquadravs_stdA[i]/mzquadravsA[i])**2+(stdmzsqsA[i]/mzsq_avsA[i])**2)
    
for i in range(lengthB):
    ULBx[i] = 1 - (mxquadravsB[i]/(3*mxsq_avsB[i]**2))
    ULBy[i] = 1 - (myquadravsB[i]/(3*mysq_avsB[i]**2))
    ULBz[i] = 1 - (mzquadravsB[i]/(3*mzsq_avsB[i]**2))
    #
    ULBdelta_x[i] = ULBx[i]*sqrt((mxquadravs_stdB[i]/mxquadravsB[i])**2+(stdmxsqsB[i]/mxsq_avsB[i])**2)
    ULBdelta_y[i] = ULBy[i]*sqrt((myquadravs_stdB[i]/myquadravsB[i])**2+(stdmysqsB[i]/mysq_avsB[i])**2)
    ULBdelta_z[i] = ULBz[i]*sqrt((mzquadravs_stdB[i]/mzquadravsB[i])**2+(stdmzsqsB[i]/mzsq_avsB[i])**2)

'''
for i in range(lengthC):
    ULCx[i] = 1 - (mxquadravsC[i]/(3*mxsq_avsC[i]**2))
    ULCy[i] = 1 - (myquadravsC[i]/(3*mysq_avsC[i]**2))
    ULCz[i] = 1 - (mzquadravsC[i]/(3*mzsq_avsC[i]**2))
    #
    ULCdelta_x[i] = ULCx[i]*sqrt((mxquadravs_stdC[i]/mxquadravsC[i])**2+(stdmxsqsC[i]/mxsq_avsC[i])**2)
    ULCdelta_y[i] = ULCy[i]*sqrt((myquadravs_stdC[i]/myquadravsC[i])**2+(stdmysqsC[i]/mysq_avsC[i])**2)
    ULCdelta_z[i] = ULCz[i]*sqrt((mzquadravs_stdC[i]/mzquadravsC[i])**2+(stdmzsqsC[i]/mzsq_avsC[i])**2)
    

for i in range(lengthD):
    ULDx[i] = 1 - (mxquadravsD[i]/(3*mxsq_avsD[i]**2))
    ULDy[i] = 1 - (myquadravsD[i]/(3*mysq_avsD[i]**2))
    ULDz[i] = 1 - (mzquadravsD[i]/(3*mzsq_avsD[i]**2))
    #
    ULDdelta_x[i] = ULDx[i]*sqrt((mxquadravs_stdD[i]/mxquadravsD[i])**2+(stdmxsqsD[i]/mxsq_avsD[i])**2)
    ULDdelta_y[i] = ULDy[i]*sqrt((myquadravs_stdD[i]/myquadravsD[i])**2+(stdmysqsD[i]/mysq_avsD[i])**2)
    ULDdelta_z[i] = ULDz[i]*sqrt((mzquadravs_stdD[i]/mzquadravsD[i])**2+(stdmzsqsD[i]/mzsq_avsD[i])**2)
    

for i in range(lengthE):
    ULEx[i] = 1 - (mxquadravsE[i]/(3*mxsq_avsE[i]**2))
    ULEy[i] = 1 - (myquadravsE[i]/(3*mysq_avsE[i]**2))
    ULEz[i] = 1 - (mzquadravsE[i]/(3*mzsq_avsE[i]**2))
    #
    ULEdelta_x[i] = ULEx[i]*sqrt((mxquadravs_stdE[i]/mxquadravsE[i])**2+(stdmxsqsE[i]/mxsq_avsE[i])**2)
    ULEdelta_y[i] = ULEy[i]*sqrt((myquadravs_stdE[i]/myquadravsE[i])**2+(stdmysqsE[i]/mysq_avsE[i])**2)
    ULEdelta_z[i] = ULEz[i]*sqrt((mzquadravs_stdE[i]/mzquadravsE[i])**2+(stdmzsqsE[i]/mzsq_avsE[i])**2)
    

for i in range(lengthF):
    ULFx[i] = 1 - (mxquadravsF[i]/(3*mxsq_avsF[i]**2))
    ULFy[i] = 1 - (myquadravsF[i]/(3*mysq_avsF[i]**2))
    ULFz[i] = 1 - (mzquadravsF[i]/(3*mzsq_avsF[i]**2))
    #
    ULFdelta_x[i] = ULFx[i]*sqrt((mxquadravs_stdF[i]/mxquadravsF[i])**2+(stdmxsqsF[i]/mxsq_avsF[i])**2)
    ULFdelta_y[i] = ULFy[i]*sqrt((myquadravs_stdF[i]/myquadravsF[i])**2+(stdmysqsF[i]/mysq_avsF[i])**2)
    ULFdelta_z[i] = ULFz[i]*sqrt((mzquadravs_stdF[i]/mzquadravsF[i])**2+(stdmzsqsF[i]/mzsq_avsF[i])**2)
    

for i in range(lengthG):
    ULGx[i] = 1 - (mxquadravsG[i]/(3*mxsq_avsG[i]**2))
    ULGy[i] = 1 - (myquadravsG[i]/(3*mysq_avsG[i]**2))
    ULGz[i] = 1 - (mzquadravsG[i]/(3*mzsq_avsG[i]**2))
    #
    ULGdelta_x[i] = ULGx[i]*sqrt((mxquadravs_stdG[i]/mxquadravsG[i])**2+(stdmxsqsG[i]/mxsq_avsG[i])**2)
    ULGdelta_y[i] = ULGy[i]*sqrt((myquadravs_stdG[i]/myquadravsG[i])**2+(stdmysqsG[i]/mysq_avsG[i])**2)
    ULGdelta_z[i] = ULGz[i]*sqrt((mzquadravs_stdG[i]/mzquadravsG[i])**2+(stdmzsqsG[i]/mzsq_avsG[i])**2)
'''


# m_x
figure()
plot(betasA, ULAx, label='L=%i'%LA)
hold('on')
plot(betasB, ULBx, label='L=%i'%LB)
#plot(betasC, ULCx, label='L=%i'%LC)
#plot(betasD, ULDx, label='L=%i'%LD)
#plot(betasE, ULEx, label='L=%i'%LE)
#plot(betasF, ULFx, label='L=%i'%LF)
#plot(betasG, ULGx, label='L=%i'%LG)
title(r'The Binder cumulant vs temperature for different L', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$U_{L,x}$', fontsize=20)
legend(loc="lower right")
show()


# m_y
figure()
plot(betasA, ULAy, label='L=%i'%LA)
hold('on')
plot(betasB, ULBy, label='L=%i'%LB)
#plot(betasC, ULCy, label='L=%i'%LC)
#plot(betasD, ULDy, label='L=%i'%LD)
#plot(betasE, ULEy, label='L=%i'%LE)
#plot(betasF, ULFy, label='L=%i'%LF)
#plot(betasG, ULGy, label='L=%i'%LG)
title(r'The Binder cumulant vs temperature for different L', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$U_{L,y}$', fontsize=20)
legend(loc="lower right")
show()

# m_z

figure()
plot(betasA, ULAz, label='L=%i'%LA)
hold('on')
plot(betasB, ULBz, label='L=%i'%LB)
#plot(betasC, ULCz, label='L=%i'%LC)
#plot(betasD, ULDz, label='L=%i'%LD)
#plot(betasE, ULEz, label='L=%i'%LE)
#plot(betasF, ULFz, label='L=%i'%LF)
#plot(betasG, ULGz, label='L=%i'%LG)
title(r'The Binder cumulant vs temperature for different L', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
legend(loc="lower right")
show()

# m_x
figure()
errorbar(betasA, ULAx, yerr=ULAdelta_x, capsize=2, label='L=%i'%LA)
hold('on')
errorbar(betasB, ULBx, yerr=ULBdelta_x, capsize=2, label='L=%i'%LB)
#errorbar(betasC, ULCx, yerr=ULCdelta_x, capsize=2, label='L=%i'%LC)
#errorbar(betasD, ULDx, yerr=ULDdelta_x, capsize=2, label='L=%i'%LD)
#errorbar(betasE, ULEx, yerr=ULEdelta_x, capsize=2, label='L=%i'%LE)
#errorbar(betasF, ULFx, yerr=ULFdelta_x, capsize=2, label='L=%i'%LF)
#errorbar(betasG, ULGx, yerr=ULGdelta_x, capsize=2, label='L=%i'%LG)
title(r'The Binder cumulant vs temperature for different L', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$U_{L,x}$', fontsize=20)
legend(loc="lower right")
show()

# m_y
figure()
errorbar(betasA, ULAy, yerr=ULAdelta_y, capsize=2, label='L=%i'%LA)
hold('on')
errorbar(betasB, ULBy, yerr=ULBdelta_y, capsize=2, label='L=%i'%LB)
#errorbar(betasC, ULCy, yerr=ULCdelta_y, capsize=2, label='L=%i'%LC)
#errorbar(betasD, ULDy, yerr=ULDdelta_y, capsize=2, label='L=%i'%LD)
#errorbar(betasE, ULEy, yerr=ULEdelta_y, capsize=2, label='L=%i'%LE)
#errorbar(betasF, ULFy, yerr=ULFdelta_y, capsize=2, label='L=%i'%LF)
#errorbar(betasG, ULGy, yerr=ULGdelta_y, capsize=2, label='L=%i'%LG)
title(r'The Binder cumulant vs temperature for different L', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$U_{L,y}$', fontsize=20)
legend(loc="lower right")
show()

# m_z
figure()
errorbar(betasA, ULAz, yerr=ULAdelta_z, capsize=2, label='L=%i'%LA)
hold('on')
errorbar(betasB, ULBz, yerr=ULBdelta_z, capsize=2, label='L=%i'%LB)
#errorbar(betasC, ULCz, yerr=ULCdelta_z, capsize=2, label='L=%i'%LC)
#errorbar(betasD, ULDz, yerr=ULDdelta_z, capsize=2, label='L=%i'%LD)
#errorbar(betasE, ULEz, yerr=ULEdelta_z, capsize=2, label='L=%i'%LE)
#errorbar(betasF, ULFz, yerr=ULFdelta_z, capsize=2, label='L=%i'%LF)
#errorbar(betasG, ULGz, yerr=ULGdelta_z, capsize=2, label='L=%i'%LG)
title(r'The Binder cumulant vs temperature for different L', fontsize=16)
xlabel(r'$\beta$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
legend(loc="lower right")
show()



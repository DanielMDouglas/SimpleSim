import os
import sys
from eventRecord import *

os.system("python3 drift.py -g rock")
print('Drifted')

evt_record = np.load('driftHits.npy', allow_pickle=True)[0]
if evt_record.inside == True:
    os.system("python3 pixelate.py")
    print('Pixelated')
    os.system("python3 reco.py")
    print('Reconstructed')

    os.system("rm driftHits.npy pixels.npy") #Delete intermediate files 
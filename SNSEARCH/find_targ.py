#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *

image = sys.argv[1]
xx = float(sys.argv[2])
yy = float(sys.argv[3])

data= loadtxt(image+'.stars')
offsets = sqrt((xx-data[:,0])**2.+(yy-data[:,1])**2.)
idx = where(offsets == min(offsets))

print 'offset = %.1f ' %offsets[idx]
print 'R_e = %2.2f, Mag = %2.2f+/-%.2f' %(data[idx][:,2], data[idx][:,3], data[idx][:,-1])


    
from msmbuilder import io
from scipy.io import *
from msmbuilder import Trajectory
import numpy
import pylab
from msmbuilder import msm_analysis

T=mmread('./l6/tProb.mtx')
map=numpy.loadtxt('./l6/Mapping.dat')
com_dist=numpy.loadtxt('Gens_com_dist.dat', usecols=(1,))
prmsd=numpy.loadtxt('Gens_p53_rmsd.dat', usecols=(1,))

frames=numpy.where(map!=-1)[0]
stateprmsd=prmsd[frames]
statecom=com_dist[frames]

eigs_m=msm_analysis.get_eigenvectors(T, 10)

#import pdb
#pdb.set_trace()

order=numpy.argsort(stateprmsd)
ordercom=numpy.argsort(statecom)

cm=pylab.cm.get_cmap('RdYlBu_r') #blue will be negative components, red positive

print numpy.shape(eigs_m[1][:,1])
print len(stateprmsd)
print len(frames)
for i in range(0,4):
    #pylab.scatter(statermsd[order], statehelix[order]*100., c=eigs_m[1][:,i], cmap=cm, s=1000*abs(eigs_m[1][:,i]), alpha=0.7)
    pylab.scatter(statecom[ordercom], stateprmsd[ordercom]*10., c=eigs_m[1][ordercom,i], cmap=cm, s=1000*abs(eigs_m[1][ordercom,i]), alpha=0.5)
    print statecom[ordercom][numpy.argmax(eigs_m[1][ordercom,i])]
    print eigs_m[1][ordercom,i][1]
#    pylab.scatter(statermsd[order], statehelix[order]*100., c=eigs_m[1][:,i], cmap=cm, s=50, alpha=0.7)
    pylab.subplots_adjust(left = 0.1, right = 1.02, bottom = 0.10, top = 0.85, wspace = 0, hspace = 0)
    CB=pylab.colorbar()
    l,b,w,h=pylab.gca().get_position().bounds
    ll, bb, ww, hh=CB.ax.get_position().bounds
    CB.ax.set_position([ll, b+0.1*h, ww, h*0.8])
    ylabel=pylab.ylabel('p53 RMSD to Bound Conformation ($\AA$)')
    xlabel=pylab.xlabel(r'p53 to S100B($\beta$$\beta$) CoM Separation ($\AA$)')
    pylab.ylim(0, max(stateprmsd)*10)
 #   loc, labels=pylab.yticks(range(0, 100, 10), [j for j in numpy.arange(0,100,10)])
#    pylab.title('Folding and Binding \n Colored by Magnitudes of Slowest Eigenvector Components')
    pylab.savefig('2deigs%i_com_prmsd.pdf' %(i),dpi=600)
    pylab.clf()
#    pylab.show()

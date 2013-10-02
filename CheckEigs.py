#!/bin/python
from msmbuilder import Trajectory, Project, io, msm_analysis, tpt
import glob
import random
from scipy.io import *
from msmbuilder import Conformation
from msmbuilder import MSMLib, msm_analysis
import optparse
import numpy
import os
from numpy import linalg
import pylab

def map_size(x):
    if x==0.5:
        size=500
    elif abs(x-0.5) < 0.1 and abs(x-0.5) > 0:
        size=300
    elif abs(x-0.5) < 0.2 and abs(x-0.5) > 0.1:
        size=200
    elif abs(x-0.5) < 0.3 and abs(x-0.5) > 0.2:
        size=100
    elif x==0:
        size=10
    elif x==1:
        size=10
    else:
        size=50
    return size

def main(dir, coarse , lag, type):
    data=dict()
    data['rmsd']=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.selfrmsd.dat' % (dir, coarse))
    #data['rmsd']=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.rmsd.dat' % (dir, coarse))
    #data['rmsd']=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.helixrmsd.dat' % (dir, coarse))
    com=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.vmd_com.dat' % (dir, coarse), usecols=(1,))
    #com=[i/com[0] for i in com]
    data['com']=com[1:]
    modeldir='%s/msml%s_coarse_r10_d%s/' % (dir, lag, coarse)
    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)

    reference=dict()
    ref_coms=numpy.loadtxt('/home/mlawrenz/p53/s100b_md/msm_80ns_p53CA/5A_cutoff/vmd_com.dat', usecols=(1,))
    ref_rmsd=numpy.loadtxt('/home/mlawrenz/p53/s100b_md/msm_80ns_p53CA/5A_cutoff/trajrmsd.dat', usecols=(1,))
    reference['rmsd']=[]
    reference['com']=[]
    for conf in range(1,10):
        n=conf-1
        reference['rmsd'].append(float(ref_rmsd[n]))
        reference['com'].append(float(ref_coms[n]))
    map_rmsd=[]
    map_com=[]
    for x in range(0, len(data['rmsd'])):
        if map[x]!=-1:
            map_com.append(data['com'][x])
            map_rmsd.append(data['rmsd'][x])
    
    map_com=numpy.array(map_com)
    map_rmsd=numpy.array(map_rmsd)
    T=mmread('%s/tProb.mtx' % modeldir)
    eigs_m=msm_analysis.get_eigenvectors(T, 10)

    order=numpy.argsort(map_rmsd)
    ordercom=numpy.argsort(map_com)

    cm=pylab.cm.get_cmap('RdYlBu_r') #blue will be negative components, red positive

    print numpy.shape(eigs_m[1][:,1])
    for i in range(0,6):
        if i==0:
            order=numpy.argsort(eigs_m[1][:,i])
            maxes=[]
            values=[]
            for n in order[::-1][:5]:
                maxes.append(numpy.where(map==n)[0])
                values.append(eigs_m[1][n,i])
            print "maxes at ",  maxes, values
        else:
            order=numpy.argsort(eigs_m[1][:,i])
            maxes=[]
            values=[]
            for n in order[::-1][:5]:
                maxes.append(numpy.where(map==n)[0])
                values.append(eigs_m[1][n,i])
            print "maxes at ",  maxes, values
            order=numpy.argsort(eigs_m[1][:,i])
            mins=[]
            values=[]
            for n in order[:5]:
                mins.append(numpy.where(map==n)[0])
                values.append(eigs_m[1][n,i])
            print "mins at ",  mins, values
        pylab.scatter(map_com[ordercom], map_rmsd[ordercom], c=eigs_m[1][ordercom,i], cmap=cm, s=1000*abs(eigs_m[1][ordercom,i]), alpha=0.5)
        print map_com[ordercom][numpy.argmax(eigs_m[1][ordercom,i])]
        print eigs_m[1][ordercom,i][1]
#       pylab.scatter(map_rmsd[order], statehelix[order]*100., c=eigs_m[1][:,i], cmap=cm, s=50, alpha=0.7)
        rdivide=max(reference['rmsd'])/10.0
        cdivide=max(reference['com'])/10.0
        pylab.plot([max(reference['com'])]*11, numpy.arange(0, max(reference['rmsd'])+rdivide/2.0, rdivide), 'r--', label='NMR Bound Ensemble')
        pylab.plot(numpy.arange(0, max(reference['com'])+cdivide, cdivide), [max(reference['rmsd'])]*11, 'r--')
        pylab.subplots_adjust(left = 0.1, right = 1.02, bottom = 0.10, top = 0.85, wspace = 0, hspace = 0)
        CB=pylab.colorbar()
        l,b,w,h=pylab.gca().get_position().bounds
        ll, bb, ww, hh=CB.ax.get_position().bounds
        CB.ax.set_position([ll, b+0.1*h, ww, h*0.8])
        CB.set_label('Eig%s Magnitudes' % i)
        ylabel=pylab.ylabel('p53 RMSD to Bound Conformation ($\AA$)')
        xlabel=pylab.xlabel(r'p53 to S100$\beta$$\beta$ Active Site COM Distance ($\AA$)')
        #pylab.ylim(0, 12)
        #pylab.xlim(0,60)
        #pylab.title('Folding and Binding \n Colored by Magnitudes of Slowest Eigenvector Components')
        pylab.legend(loc=8, frameon=False)
        pylab.savefig('%s/2deigs%i_com_prmsd.png' %(modeldir, i),dpi=300)
        pylab.show()
        #pylab.clf()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    parser.add_option('-c', '--coarse', dest='coarse',
                      help='coarse grain cutoff')
    parser.add_option('-l', '--lag', dest='lag',
                      help='lag time')
    parser.add_option('-t', '--type', dest='type',
                      help='type')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(dir=options.dir, coarse=options.coarse, lag=options.lag, type=options.type)

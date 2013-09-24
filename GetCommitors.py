#!/bin/python
from msmbuilder import Trajectory, Project, io, msm_analysis, tpt
import glob
import random
from scipy.io import *
from msmbuilder import Conformation
from msmbuilder import MSMLib
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
    rmsd=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.rmsd.dat' % (dir, coarse), usecols=(2,))
    data['rmsd']=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.selfrmsd.dat' % (dir, coarse))
    com=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.vmd_com.dat' % (dir, coarse), usecols=(1,))
    com=[i/com[0] for i in com]
    data['com']=com[1:]
    modeldir='%s/msml%s_coarse_r10_d%s/' % (dir, lag, coarse)
    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)

    map_rmsd=[]
    map_com=[]
    for x in range(0, len(data['rmsd'])):
        if map[x]!=-1:
            map_com.append(data['com'][x])
            map_rmsd.append(data['rmsd'][x])
    
    map_com=numpy.array(map_com)
    map_rmsd=numpy.array(map_rmsd)
    T=mmread('%s/tProb.mtx' % modeldir)
    unbound=numpy.loadtxt('%s/tpt-%s/unbound_%s_states.txt' % (modeldir, type, type), dtype=int)
    bound=numpy.loadtxt('%s/tpt-%s/bound_%s_states.txt' % (modeldir, type, type), dtype=int)

    Tdense=T.todense()
    data=dict()
    for i in unbound:
        for j in unbound:
            if Tdense[i,j]!=0:
                if i not in data.keys():
                    data[i]=[]
                data[i].append(j)
    print data
    cm=pylab.cm.get_cmap('RdYlBu_r') #blue will be negative components, red positive
    Q=tpt.calculate_committors(unbound, bound, T)
    for i in range(0,len(Q)):
        if Q[i]>0.40 and Q[i]<0.6:
            print i
    pylab.scatter(map_com, map_rmsd, c=Q, cmap=cm, alpha=0.7, s=[map_size(i) for i in Q])
    pylab.colorbar()
    pylab.show()


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

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
    else:
        size=50
    return size

def main(modeldir,  type):
    data=dict()
    data['rmsd']=numpy.loadtxt('%s/Gens.rmsd.dat' % modeldir, usecols=(2,))
    com=numpy.loadtxt('%s/Gens.vmd_com.dat' % modeldir, usecols=(1,))
    refcom=com[0]
    data['com']=com[1:]
    data['com']=numpy.array(data['com'])
    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)
    frames=numpy.where(map!=-1)[0]

    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)
    unbound=numpy.loadtxt('%s/tpt-%s/unbound_%s_states.txt' % (modeldir, type, type), dtype=int)
    bound=numpy.loadtxt('%s/tpt-%s/bound_%s_states.txt' % (modeldir, type, type), dtype=int)

    paths=io.loadh('%s/tpt-%s/Paths.h5' % (modeldir, type))

    committors=numpy.loadtxt('%s/commitor_states.txt' % modeldir, dtype=int)
    print committors
    colors=['red', 'orange', 'green', 'cyan', 'blue', 'purple']
    colors=colors*40
    map_com=[]
    map_rmsd=[]
    for x in range(0, len(data['rmsd'])):
        if map[x]!=-1:
            map_com.append(data['com'][x])
            map_rmsd.append(data['rmsd'][x])
    map_com=numpy.array(map_com)
    map_rmsd=numpy.array(map_rmsd)
    if type=='strict':
        ref=5
    elif type=='medium':
        ref=10
    elif type=='loose':
        ref=15
    for p in range(0, 10):
        pylab.figure()
        path=paths['Paths'][p]
        print "Bottleneck", paths['Bottlenecks'][p]
        flux=paths['fluxes'][p]/paths['fluxes'][0]
        if flux < 0.2:
            break
        print "flux %s" % flux
        frames=numpy.where(path!=-1)[0]
        path=numpy.array(path[frames], dtype=int)
        size=(paths['fluxes'][p]/paths['fluxes'][0])*1000
        for j in paths['Bottlenecks'][p]:
            pylab.scatter(map_com[j], map_rmsd[j], marker='x', c=colors[p], alpha=0.7, s=size*2)
            location=numpy.where(committors==j)[0]
            if location.size:
                print "path %s state %s bottleneck in committors" % (p, j)
                print map_com[j], map_rmsd[j]
        pylab.scatter(map_com[path], map_rmsd[path], c=colors[p], alpha=0.7, s=size)
        pylab.plot([ref]*10, numpy.arange(0, max(data['rmsd']), max(data['rmsd'])/10), 'r--')
        pylab.hold(True)
        pylab.xlabel('PL-COM')
        pylab.ylabel('L RMSD')
        pylab.xlim(0,max(map_com)+5)
        pylab.ylim(0,max(map_rmsd+5))
    pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    parser.add_option('-t', '--type', dest='type',
                      help='type')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(modeldir=options.dir, type=options.type)


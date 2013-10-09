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
    else:
        size=50
    return size

def main(modeldir, type):
    data=dict()
    project=Project.load_from('%s/ProjectInfo.yaml' % modeldir.split('Data')[0])
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

    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)
    T=mmread('%s/tProb.mtx' % modeldir)
    paths=io.loadh('%s/tpt-%s/Paths.h5' % (modeldir, type))

    steps=10000
    for startstate in unbound:
        print "on start state %s" % startstate
        movie=project.empty_traj()
        traj=msm_analysis.sample(T, int(startstate),int(steps))
        for (n, state) in enumerate(traj):
            t=project.get_random_confs_from_states(ass['arr_0'], [int(state),], 20)
            if n==0:
                movie['XYZList']=t[0]['XYZList']
            else:
                movie['XYZList']=numpy.vstack((movie['XYZList'], t[0]['XYZList']))
        movie.save_to_xtc('%s/tpt-%s/movie_state%s_100microsec.xtc' % (modeldir, type, startstate))

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
    main(dir=options.dir,  type=options.type)


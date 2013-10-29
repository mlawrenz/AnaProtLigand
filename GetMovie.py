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

def main(modeldir, start, type):
    start=int(start)
    data=dict()
    project=Project.load_from('%s/ProjectInfo.yaml' % modeldir.split('Data')[0])
    files=glob.glob('%s/fkbp*xtal.pdb' % modeldir.split('Data')[0])
    pdb=files[0]
    unbound=numpy.loadtxt('%s/tpt-%s/unbound_%s_states.txt' % (modeldir, type, type), dtype=int)
    T=mmread('%s/tProb.mtx' % modeldir)
    startstate=unbound[start]
    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)

    steps=100000
    print "on start state %s" % startstate
    if os.path.exists('%s/tpt-%s/movie_state%s_1millisec.states.dat' % (modeldir, type, startstate)):
        print "loading from states"
        traj=numpy.loadtxt('%s/tpt-%s/movie_state%s_1millisec.states.dat' % (modeldir, type, startstate))
    else:
        traj=msm_analysis.sample(T, int(startstate),int(steps))
        numpy.savetxt('%s/tpt-%s/movie_state%s_1millisec.states.dat' % (modeldir, type, startstate), traj)
    print "checking for chkpt file"
    checkfile=glob.glob('%s/tpt-%s/movie_state%s_*chkpt' % (modeldir, type, startstate))
    if len(checkfile) > 0:
        movie=Trajectory.load_from_xtc(checkfile[0], PDBFilename=pdb)
        n=int(checkfile[0].split('xtc.state')[1].split('chkpt')[0])
        os.system('mv %s %s.chkpt.cp' % (checkfile[0], checkfile[0].split('.xtc')[0]))
        print "checkpointing at state index %s out of %s" % (n, len(traj))
        checkfile=checkfile[0]
        restart=True
    else:
        restart=False
        n=0
        movie=project.empty_traj()
    while n < len(traj):
        print "on state %s" % n
        state=int(traj[n])
        t=project.get_random_confs_from_states(ass['arr_0'], [int(state),], 10)
        if n==0:
            movie['XYZList']=t[0]['XYZList']
            n+=1
            continue
        elif n % 100==0:
            movie['XYZList']=numpy.vstack((movie['XYZList'], t[0]['XYZList']))
            if restart==True:
                os.system('mv %s %s.chkpt.cp' % (checkfile, checkfile.split('.xtc')[0]))
            movie.save_to_xtc('%s/tpt-%s/movie_state%s_1millisec.xtc.state%schkpt' % (modeldir, type, startstate, n))
            checkfile='%s/tpt-%s/movie_state%s_1millisec.xtc.state%schkpt' % (modeldir, type, startstate, n)
            n+=1
            continue
        elif n!=0:
            movie['XYZList']=numpy.vstack((movie['XYZList'], t[0]['XYZList']))
            n+=1
            continue
    movie.save_to_xtc('%s/tpt-%s/movie_state%s_1millisec.xtc' % (modeldir, type, startstate))

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    parser.add_option('-s', '--start', dest='start',
                      help='start index of unbound states')
    parser.add_option('-t', '--type', dest='type',
                      help='type')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(modeldir=options.dir,  start=options.start, type=options.type)


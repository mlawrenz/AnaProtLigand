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

def main(dir, coarse , lag, type):
    data=dict()
    modeldir='%s/msml%s_coarse_r10_d%s/' % (dir, lag, coarse)

    project = Project.load_from('ProjectInfo.yaml')
    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)
    T=mmread('%s/tProb.mtx' % modeldir)
    paths=io.loadh('%s/tpt-%s/Paths.h5' % (modeldir, type))

    if not os.path.exists('%s/adaptive-states/' % modeldir):
        os.mkdir('%s/adaptive-states/' % modeldir)
    for state in sorted(set(ass['arr_0'].flatten())):
        if state!=-1:
            t=project.get_random_confs_from_states(ass['arr_0'], [int(state),], 5)
            for i in range(0, 5):
                print state, i
                (a, b, c) =t[0]['XYZList'].shape
                movie=project.empty_traj()
                movie['XYZList']=numpy.zeros((1, b, c), dtype=numpy.float32)
                movie['XYZList'][0]=t[0]['XYZList'][i]
                movie.save_to_pdb('%s/adaptive-states/state%s-%s.pdb' % (modeldir, int(state), i))

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


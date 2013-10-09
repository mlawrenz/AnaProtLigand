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

def main(modeldir):
    data=dict()
    project=Project.load_from('%s/ProjectInfo.yaml' % modeldir.split('Data')[0])
    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)
    T=mmread('%s/tProb.mtx' % modeldir)

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
    parser.add_option('-d', '--modeldir', dest='modeldir',
                      help='model directory')
    parser.add_option('-t', '--type', dest='type',
                      help='type')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(modeldir=options.modeldir, type=options.type)


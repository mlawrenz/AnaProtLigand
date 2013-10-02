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
    project = Project.load_from('ProjectInfo.yaml')
    modeldir='%s/msml%s_coarse_r10_d%s/' % (dir, lag, coarse)
    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)
    T=mmread('%s/tProb.mtx' % modeldir)
    paths=io.loadh('%s/tpt-%s/Paths.h5' % (modeldir, type))
    
    for p in range(0, 20):
        movie=project.empty_traj()
        path=paths['Paths'][p]
        flux=paths['fluxes'][p]/paths['fluxes'][0]
        if flux < 0.2:
            break
        print "flux %s" % flux
        frames=numpy.where(path!=-1)[0]
        path=numpy.array(path[frames], dtype=int)
        for (n, state) in enumerate(path):
            t=project.get_random_confs_from_states(ass['arr_0'], [int(state),], 20)
            if n==0:
                movie['XYZList']=t[0]['XYZList']
            else:
                movie['XYZList']=numpy.vstack((movie['XYZList'], t[0]['XYZList']))
        movie.save_to_xtc('%s/tpt-%s/path%s_sample20.xtc' % (modeldir, type, p))

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


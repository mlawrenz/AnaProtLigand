#!/bin/python 

import numpy
import os
import pylab
from scipy.sparse import linalg
from scipy.io import *
from msmbuilder import msm_analysis, Project, io
import optparse


def main(dir, lag):
    project = Project.load_from('ProjectInfo.yaml')
    T=mmread('%s/tProb.mtx' % dir)
    map=numpy.loadtxt('%s/Mapping.dat' % dir)
    ass=io.loadh('%s/Assignments.h5' % dir)
    frames=numpy.where(map!=-1)[0]
    trim=numpy.where(map==-1)[0]
    rmsd=rmsd[frames]
    #confs=project.get_random_confs_from_states(ass['arr_0'], frames, 5)
    #for n in range(0, len(confs)):
    #    confs[n].save_to_xtc('%s/MSMl%s/state%s.xtc' % (dir, lag, frames[n]))
    print "%s trimmed %s frames, max rmsd %s" % (dir, len(trim), max(rmsd)) 
    eigs=msm_analysis.get_eigenvectors(T, 20)
    order=numpy.argsort(rmsd)
    for eigen in eigs[0][1:]:
        implied = -float(lag)/numpy.log(eigen)
        pylab.plot(range(0,10), [numpy.log(implied/10.0)]*10)
    pylab.title('eigen spectrum %s' % dir)
    pylab.show()
    for i in range(0, 10):
        print eigs[0][i]
        print "sum is %s" % numpy.sum(eigs[1][order,i])
        if len(rmsd)!=len(eigs[1][:,i]):
            import pdb
            pdb.set_trace()
        pylab.plot(range(0, len(eigs[1][order,i])), eigs[1][order,i])
        pylab.hold(False)
        loc,labels=pylab.xticks(range(0, len(eigs[1][order,i]), 5), [("%0.2f" % x) for x in rmsd[order][::5]])
        #pylab.xlim(-20, len(eigs[1][order,i]))
        pylab.xlabel('Ligand RMSD of Eigenvector Components to Bound State')
        pylab.ylabel('Eigenvector %s' % (i+1))
        if not os.path.exists('%s/eigs' % dir):
            os.mkdir('%s/eigs' % dir)
        pylab.savefig('%s/eigs/mle-eig%i.png' % (dir, i))
        pylab.title('%s' % dir)
        pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    parser.add_option('-l', '--lag', dest='lag',
                      help='lag time')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(dir=options.dir, lag=options.lag)


from msmbuilder import io
from msmbuilder import MSMLib, Project
from scipy.io import *
from scipy.sparse import lil_matrix
import os
import scipy
import optparse
import pickle
import numpy
import pylab

def main(assfile, lag, itstart=0):
    lag=int(lag)
    Assignments=io.loadh(assfile)
    dir=os.path.dirname(assfile)
    newdir='%s/sample-counts' % dir
    proj=Project.load_from('%s/ProjectInfo.yaml' % dir.split('Data')[0])
    multinom=sum(proj.traj_lengths)
    if not os.path.exists(newdir):
        os.mkdir(newdir)
    if 'Data' in Assignments.keys():
        Assignments=Assignments['Data']
    else:
        Assignments=Assignments['arr_0']
    print Assignments.shape
    NumStates = max(Assignments.flatten()) + 1
    Counts = MSMLib.get_count_matrix_from_assignments(Assignments, lag_time=int(lag), sliding_window=True)
    Counts=Counts.todense()
    Counts=Counts*(1.0/lag)
    T=numpy.array(Counts)
    for (i,j) in zip(numpy.where(T==0)[0], numpy.where(T==0)[1]):
        T[i,j]=1    
    Popsample=dict()
    iteration=itstart
    print "iterating thru tCount samples"
    while iteration < 100:
        print "sampling iteration %s" % iteration
        newT=scipy.sparse.lil_matrix((int(NumStates),int(NumStates)),dtype='float32') 
        for i in range(0, T.shape[1]):
            transitions = numpy.row_stack((numpy.array([i]*NumStates),numpy.arange(0, NumStates)))
            pvals=numpy.array([x/sum(T[i]) for x in T[i]]) 
            counts=numpy.random.multinomial(int(multinom), pvals, size=1)
            newT=newT+scipy.sparse.coo_matrix((counts[0], transitions),shape=(NumStates,NumStates))
        rev_counts, t_matrix, Populations, Mapping = MSMLib.build_msm(newT, symmetrize='MLE', ergodic_trimming=True)
        scipy.io.mmwrite('%s/tProb-%s' % (newdir, iteration), t_matrix)
        iteration+=1
    

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-a', '--assfile', dest='assfile',
                      help='input  Assignments.h5 file')
    parser.add_option('-l', '--lag', dest='lag',
                      help='lag time')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(assfile=options.assfile, lag=options.lag)

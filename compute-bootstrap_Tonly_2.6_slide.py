from msmbuilder import io
from msmbuilder import MSMLib
from scipy.io import *
from scipy.sparse import lil_matrix
import os
import scipy
import optparse
import pickle
import numpy
import pylab

def main(dir, lag, multinom, itstart):
    itstart=int(itstart)
    lag=int(lag)
    multinom=int(multinom)
    Assignments=io.loadh('%s/Assignments.h5' % dir)
    if 'Data' in Assignments.keys():
        Assignments=Assignments['Data']
    else:
        Assignments=Assignments['arr_0']
    print Assignments.shape
    NumStates = max(Assignments.flatten()) + 1
    Counts = MSMLib.get_count_matrix_from_assignments(Assignments, lag_time=int(lag), sliding_window=True)
    Counts=Counts.todense()
    Counts=Counts*(1.0/lag)
    #T=mmread('newT.mtx')
    #T=mmread('%s/MSMl%s/tCounts.mtx' % (dir, lag))
    T=numpy.array(Counts)
    for (i,j) in zip(numpy.where(T==0)[0], numpy.where(T==0)[1]):
        T[i,j]=1    
    Popsample=dict()
    iteration=itstart
    print "iterating thru tCount samples"
    while iteration < itstart+5:
        print "sampling iteration %s" % iteration
        newT=scipy.sparse.lil_matrix((int(NumStates),int(NumStates)),dtype='float32') 
        for i in range(0, T.shape[1]):
            transitions = numpy.row_stack((numpy.array([i]*NumStates),numpy.arange(0, NumStates)))
            pvals=numpy.array([x/sum(T[i]) for x in T[i]]) 
            counts=numpy.random.multinomial(int(multinom), pvals, size=1)
            newT=newT+scipy.sparse.coo_matrix((counts[0], transitions),shape=(NumStates,NumStates))
        rev_counts, t_matrix, Populations, Mapping = MSMLib.build_msm(newT, symmetrize='MLE', ergodic_trimming=True)
        scipy.io.mmwrite('%s/tProb-%s' % (dir, iteration), t_matrix)
        iteration+=1
    

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='input dir with Assignments.h5 file')
    parser.add_option('-l', '--lag', dest='lag',
                      help='lag time')
    parser.add_option('-i', '--itmax', dest='itmax',
                      help='tCount bootstrap iterations (N for standard deviation)')
    parser.add_option('-m', '--multinom', dest='multinom',
            help='size of multinomian distribution sample: recommend # of trajectory snapshots in project')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(dir=options.dir, lag=options.lag, itstart=options.itmax, multinom=options.multinom)

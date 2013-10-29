import multiprocessing
import time
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

def parallel_get_matrix(input):
    print "working"
    (Ttest, multinom, n)=input
    NumStates=Ttest.shape[0]
    numpy.random.seed(int(time.time()*(n+1)))
    newT=scipy.sparse.lil_matrix((int(NumStates),int(NumStates)),dtype='float32')
    for i in range(0, Ttest.shape[1]):
        transitions = numpy.row_stack((numpy.array([i]*NumStates),numpy.arange(0, NumStates)))
        pvals=numpy.array([x/sum(Ttest[i]) for x in Ttest[i]])
        counts=numpy.random.multinomial(int(multinom), pvals, size=1)
        newT=newT+scipy.sparse.coo_matrix((counts[0], transitions),shape=(NumStates,NumStates))
    frames=numpy.where(newT==0)
    newT[frames]=1
    return newT

def main(assfile, lag, nproc):
    lag=int(lag)
    nproc=int(nproc)
    Assignments=io.loadh(assfile)
    num=int(assfile.split('Assignments_sub')[1].split('.h5')[0])
    dir=os.path.dirname(assfile)
    newdir='%s/boot-sub%s' % (dir, num)
    ref_sub=numpy.loadtxt('%s/times.h5' % dir, usecols=(1,))
    ref_total=numpy.loadtxt('%s/times.h5' % dir, usecols=(2,))
    times=dict()
    for (i,j) in zip(ref_sub, ref_total):
        times[i]=j

    proj=Project.load_from('%s/ProjectInfo.yaml' % dir.split('Data')[0])
    multinom=int(times[num])
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
    frames=numpy.where(T==0)
    T[frames]=1
    Popsample=dict()
    iteration=0
    total_iteration=100/nproc
    print "%s total iterations" % total_iteration
    if 100 % nproc != 0:
        remain=100 % nproc
    else:
        remain=False
    print "iterating thru tCount samples"
    count=0
    while iteration < 100:
        if count*nproc > 100:
            nproc=remain
        print "sampling iteration %s" % iteration
        Tfresh=T.copy()
        input = zip([Tfresh]*nproc, [multinom]*nproc, range(0, NumStates))
        pool = multiprocessing.Pool(processes=nproc)
        result = pool.map_async(parallel_get_matrix, input)
        result.wait()
        all = result.get()
        pool.terminate()
        for c_matrix in all:
            scipy.io.mmwrite('%s/tCounts-%s' % (newdir, iteration), c_matrix)
            #rev_counts, t_matrix, Populations, Mapping=x
            #scipy.io.mmwrite('%s/tProb-%s' % (newdir, iteration), t_matrix)
            #numpy.savetxt('%s/Populations-%s' % (newdir, iteration), Populations)
            #numpy.savetxt('%s/Mapping-%s' % (newdir, iteration), Mapping)
            iteration+=1
        count+=1
        print "dont with iteration %s" % iteration*nproc
    

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-a', '--assfile', dest='assfile',
                      help='input  Assignments.h5 file')
    parser.add_option('-l', '--lag', dest='lag',
                      help='lag time')
    parser.add_option('-n', '--nproc', dest='nproc',
                      help='num processors')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(assfile=options.assfile, lag=options.lag, nproc=options.nproc)

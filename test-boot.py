import multiprocessing
import time
import random
from scipy.sparse import lil_matrix
from msmbuilder import io
from msmbuilder import MSMLib, Project
from scipy.io import *
import os
import scipy
import optparse
import pickle
import numpy
import pylab



def parallel_get_matrix(input):
    (T, multinom, n)=input
    #print "working"
    NumStates=T.shape[0]
    numpy.random.seed(int(time.time()*(n+1)))
    newT=scipy.sparse.lil_matrix((int(NumStates),int(NumStates)),dtype='float32')
    for i in range(0, T.shape[1]):
        transitions = numpy.row_stack((numpy.array([i]*NumStates),numpy.arange(0, NumStates)))
        pvals=numpy.array([x/sum(T[i]) for x in T[i]])
        counts=numpy.random.multinomial(int(multinom), pvals, size=1)
        newT=newT+scipy.sparse.coo_matrix((counts[0], transitions),shape=(NumStates,NumStates))
    return newT.todense()

def parallel_get_matrix_multi_row(input):
    #print "working"
    (T, multinom, rows)=input
    NumStates=T.shape[0]
    tmp=scipy.sparse.lil_matrix((int(NumStates),int(NumStates)))
    for i in rows:
        numpy.random.seed(int(time.time()*(i+1)))
        NumStates=T.shape[0]
        transitions = numpy.row_stack((numpy.array([i]*NumStates),numpy.arange(0, NumStates)))
        pvals=numpy.array([x/sum(T[i]) for x in T[i]])
        counts=numpy.random.multinomial(int(multinom), pvals, size=1)
        tmp=tmp+scipy.sparse.coo_matrix((counts[0], transitions),shape=(NumStates,NumStates))
    return tmp

def parallel_get_matrix_row(input):
    #print "working"
    (T, multinom, i)=input
    numpy.random.seed(int(time.time()*(i+1)))
    NumStates=T.shape[0]
    transitions = numpy.row_stack((numpy.array([i]*NumStates),numpy.arange(0, NumStates)))
    pvals=numpy.array([x/sum(T[i]) for x in T[i]])
    counts=numpy.random.multinomial(int(multinom), pvals, size=1)
    return scipy.sparse.coo_matrix((counts[0], transitions),shape=(NumStates,NumStates))

def test_map(T, multinom):
    #look over rows
    NumStates=T.shape[0]
    newT=scipy.sparse.lil_matrix((int(NumStates),int(NumStates)),dtype='float32')
    rows=dict()
    for i in range(0, T.shape[1]):
        input=zip( [T], [multinom,], [i,])
        result = map(parallel_get_matrix_row, input)
        rows[i]=result[0]
        #newT=newT+result[0]
    for i in sorted(rows.keys()):
        newT=newT+rows[i]
    return newT.todense()




def compute_multi_rows(input, nproc):
    pool = multiprocessing.Pool(processes=nproc)
    result = pool.map_async(parallel_get_matrix_multi_row, input)
    result.wait()
    all = result.get()
    pool.terminate()
    return all

def test_pool_map(T, multinom, nproc):
    #look over rows
    NumStates=T.shape[0]
    newT=scipy.sparse.lil_matrix((int(NumStates),int(NumStates)))
    rows=[]
    remain=T.shape[0] % nproc
    row_indices=range(0, T.shape[1])
    chunks=dict()
    count=0
    chunksize=T.shape[1]/nproc
    chunks[count]=[]
    remain=T.shape[1] % nproc
    remainmax=T.shape[1]-remain
    for i in row_indices:
        if i % chunksize==0:
            if i==0:
                chunks[count].append(i)
                continue
            else:
                count+=1
                chunks[count]=[]
                chunks[count].append(i)
        else:
            chunks[count].append(i)
    ordered_chunks=[]
    for key in sorted(sorted(chunks.keys())):
        if key >= nproc:
            remainchunk=chunks[key]
            chunks.pop(key)
        else:
            ordered_chunks.append(chunks[key])
            remainchunk=None
    # distribute to processors
    input = zip([T]*nproc, [multinom]*nproc, ordered_chunks)
    all=compute_multi_rows(input, nproc)
    for (n, k) in enumerate(all):
        newT=newT+all[n]
    # do remain
    if remainchunk!=None:
        for i in remainchunk:
            input=zip( [T], [multinom,], [i,])
            result = map(parallel_get_matrix_row, input)
            newT=newT+result[0]
    else:
        pass
    return newT.todense()



def main(nproc):
    multinom=1000000
    nproc=int(nproc)
    NumStates=10
    C = scipy.sparse.lil_matrix((int(NumStates), int(NumStates)))
    C=C.todense()
    for i in range(0, NumStates):
        for j in range(0, NumStates):
            C[i,j]=random.randint(1,100)
    ### testing non parallel
    #serial=dict()
    #for iteration in range(0,4):
    #    T=numpy.array(C.copy())
    #    newT=test_map(T, multinom)
    #    serial[iteration]=newT
    ### testing all at once
    all_matrices=dict()
    start = time.clock()
    print "on all-at-once method"
    T=numpy.array(C.copy())
    input = zip([T]*nproc, [multinom]*nproc, range(0, nproc) )
    iteration=dict()
    pool = multiprocessing.Pool(processes=nproc)
    result = pool.map_async(parallel_get_matrix, input)
    result.wait()
    all = result.get()
    pool.close()
    pool.join()
    elapsed = (time.clock() - start)
    print "parallel T %s elapsed s" % elapsed
    for (n,z) in enumerate(all):
        all_matrices[n]=z
        #print "all method iteration %s last row of final matrix: " % n, z[-1]
    ## testing row by row
    print "on row-by-row method"
    start = time.clock()
    row_matrices=dict()
    for iteration in range(0,4):
        T=numpy.array(C.copy())
        newT=test_pool_map(T, multinom, nproc)
        frames=numpy.where(newT==0)
        newT[frames]=1
        #print "row method iteration %s last row of final matrix: " % iteration, newT[-1]
        #print newT
        #print "-----------------"
        #row_matrices[iteration]=newT
    elapsed = (time.clock() - start)
    print "parallel_row %s elapsed s" % elapsed
    #print "old T", C
    #print "new T", newT
    

if __name__ == "__main__":
    main(4)


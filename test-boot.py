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

global counter
counter = 0

def cb(r):
    global counter
    print counter, r
    counter +=1
    
def parallel_get_matrix(input):
    print "working"
    (Ttest, multinom)=input
    NumStates=Ttest.shape[0]
    newT=scipy.sparse.lil_matrix((int(NumStates),int(NumStates)),dtype='float32')
    for i in range(0, Ttest.shape[1]):
        transitions = numpy.row_stack((numpy.array([i]*NumStates),numpy.arange(0, NumStates)))
        pvals=numpy.array([x/sum(Ttest[i]) for x in Ttest[i]])
        counts=numpy.random.multinomial(int(multinom), pvals, size=1)
        newT=newT+scipy.sparse.coo_matrix((counts[0], transitions),shape=(NumStates,NumStates))
    return newT.todense()

def parallel_get_matrix_row(input):
    print "working"
    (T, multinom, i)=input
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



def compute_rows_all(input, nproc):
    pool = multiprocessing.Pool(processes=nproc)
    result = pool.map_async(parallel_get_matrix, input)
    result.wait()
    all = result.get()
    pool.terminate()
    return all

def compute_rows(input, nproc):
    pool = multiprocessing.Pool(processes=nproc)
    result = pool.map_async(parallel_get_matrix_row, input)
    result.wait()
    all = result.get()
    pool.terminate()
    return all

def test_pool_map(T, multinom, nproc):
    #look over rows
    NumStates=T.shape[0]
    newT=scipy.sparse.lil_matrix((int(NumStates),int(NumStates)),dtype='float32')
    rows=[]
    remain=T.shape[0] % nproc
    row_indices=range(0, T.shape[1])
    for i in row_indices:
        if i==0:
            rows.append(i)
            continue
        # below is if at end without remainder
        if i==max(row_indices):
            rows.append(i)
            input = zip([T]*nproc, [multinom]*nproc, rows)
            all=compute_rows(input, nproc)
            for (n, k) in enumerate(rows):
                newT=newT+all[n]
            break
        # below is if at end with remainder
        if i == T.shape[0]-remain:
            if len(rows)>0:
                # finish previous list triggered
                input = zip([T]*nproc, [multinom]*nproc, rows)
                all=compute_rows(input, nproc)
                for (n, k) in enumerate(rows):
                    newT=newT+all[n]
                    rows=[]
            nproc=remain
            for k in range(i, i+remain):
                rows.append(k)
            input = zip([T]*nproc, [multinom]*nproc, rows)
            all=compute_rows(input, nproc)
            for (n, k) in enumerate(rows):
                newT=newT+all[n]
            break
        if i % nproc == 0:
            input = zip([T]*nproc, [multinom]*nproc, rows)
            all=compute_rows(input, nproc)
            for (n, k) in enumerate(rows):
                newT=newT+all[n]
            rows=[]
            rows.append(i)
        else:
            rows.append(i)
    return newT.todense()



def main(nproc):
    multinom=100
    nproc=int(nproc)
    NumStates=15
    C = scipy.sparse.lil_matrix((int(NumStates), int(NumStates)))
    C=C.todense()
    for i in range(0, NumStates):
        for j in range(0, NumStates):
            C[i,j]=random.randint(1,100)
    ### testing non parallel
    serial=dict()
    for iteration in range(0,4):
        T=numpy.array(C.copy())
        newT=test_map(T, multinom)
        serial[iteration]=newT
    import pdb
    pdb.set_trace()
    ### testing all at once
    all_matrices=dict()
    start = time.clock()
    T=numpy.array(C.copy())
    input = zip([T]*nproc, [multinom]*nproc)
    all=compute_rows_all(input, nproc)
    for (n,z) in enumerate(all):
        all_matrices[n]=z
    elapsed = (time.clock() - start)
    print "parallel T %s elapsed s" % elapsed
    ### testing row by row
    row_matrices=dict()
    for iteration in range(0,4):
        T=numpy.array(C.copy())
        start = time.clock()
        newT=test_pool_map(T, multinom, nproc)
        row_matrices[iteration]=newT
    elapsed = (time.clock() - start)
    print "parallel_row %s elapsed s" % elapsed
    import pdb
    pdb.set_trace()
    #print "old T", C
    #print "new T", newT
    

if __name__ == "__main__":
    main(4)


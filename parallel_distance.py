import multiprocessing
import pickle
from numpy import *

def parallel_get_distance(input):
    (ligand_coors, protein_coors)=input
    minval=100000
    for ligand_coor in ligand_coors:
        for coor in protein_coors:
            val=sqrt(((ligand_coor[0]-coor[0])**2+(ligand_coor[1]-coor[1])**2+(ligand_coor[2]-coor[2])**2))
            if val < minval:
                minval=val
    return minval


def parallel_distance(nproc, ligcoors, protcoors, ofile, completed=None):
    nproc=int(nproc)
    maxstates=len(set(ligcoors.keys()))
    if completed==None:
        completed=dict()
    count=0
    ligand_coor_chunks=[]
    protein_coor_chunks=[]
    need_states=[]
    targets=[x for x in range(0, maxstates) if x not in completed.keys()]
    if maxstates % nproc != 0:
        remain=maxstates % nproc
        remain=targets[::-1][:remain]
    else:
        remain=[-100,]
    r=False
    for state in [x for x in range(0, maxstates) if x not in completed.keys()]:
        if state not in completed.keys():
            ligand_coor_chunks.append(ligcoors[state])
            protein_coor_chunks.append(protcoors[state])
            need_states.append(state)
            count+=1
            if state >= min(remain):
                r=True
                count=nproc
        if count==nproc:
            if r==True:
                nproc=len(remain)
            count=0
            input = zip(ligand_coor_chunks, protein_coor_chunks)
            pool = multiprocessing.Pool(processes=nproc)
            result = pool.map_async(parallel_get_distance, input)
            result.wait()
            minvals = result.get()
            pool.terminate()
            for (n, state) in enumerate(need_states):
                print "completed state %s" % state
                completed[state]=minvals[n]
            ohandle=open(ofile, 'wb')
            pickle.dump(completed, ohandle)
            ohandle.close()
            need_states=[]
            ligand_coor_chunks=[]
            protein_coor_chunks=[]
    return completed



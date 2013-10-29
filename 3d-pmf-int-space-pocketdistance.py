from msmbuilder import Trajectory, PMF3D
import multiprocessing
import optparse
import pylab
from numpy import * 
import glob
import os
import sys
import pickle
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
from scipy import integrate

def get_minmax(ligcoors, map):
    # build grid of min/max ligand coords from Gens.vmd_ligcoords.dat
    mapped_ligcoors=dict()
    n=0
    for i in sorted(ligcoors.keys()):
        if map[i]!=-1:
            mapped_ligcoors[n]=ligcoors[i]
            n+=1
    xmin=100000
    xmax=0
    ymin=100000
    ymax=0
    zmin=100000
    zmax=0
    dx=1
    dy=1
    dz=1
    coordinates=['x', 'y', 'z']
    mins=[xmin, ymin, zmin]
    maxes=[xmax, ymax, zmax]
    for i in sorted(mapped_ligcoors.keys()):
        for j in range(0, len(mapped_ligcoors[i])):
            for (n, k) in enumerate(coordinates):
                if mapped_ligcoors[i][j][n]<=mins[n]:
                    mins[n]=mapped_ligcoors[i][j][n]
                elif mapped_ligcoors[i][j][n]>=maxes[n]:
                    maxes[n]=mapped_ligcoors[i][j][n]
    lengths=dict()
    for (n, i) in enumerate(coordinates):
        lengths[n]=int(((round(maxes[n])+1)-round(mins[n]))/1.0)
    box_volume=lengths[0]*dx*lengths[1]*dy*lengths[2]*dz
    print "actual box volume %s" % box_volume
    print "input box volume %s" % ((6.79133*10)*(6.79133*10)*(4.80219*10))
    total=max(lengths.values())
    ranges=dict()
    for (n, i) in enumerate(coordinates):
        pad=-(lengths[n]-total)
        ranges[n]=range(int(round(mins[n])), int(round(maxes[n]+1+(pad*1.0+1))), 1)
    return mapped_ligcoors, ranges[0], ranges[1], ranges[2], box_volume


def parallel_distance(nproc, ligcoors, protcoors, ofile, completed=None):
    nproc=int(nproc)
    maxstates=len(set(ligcoors.keys()))
    if completed==None:
        completed=dict()
    count=0
    ligand_coor_chunks=[]
    protein_coor_chunks=[]
    need_states=[]
    if maxstates % nproc != 0:
        remain=maxstates % nproc
    targets=[x for x in range(0, maxstates) if x not in completed.keys()]
    remain=targets[::-1][:remain]
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
            result = pool.map_async(PMF3D.parallel_get_distance, input)
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

def main(system, genfile, lag, volume, nproc):
    lag=int(lag)
    dir=os.path.dirname(genfile)
    if not os.path.exists('%s/target.txt' % dir.split('Data')[0]):
        print "need target.txt with experimental value"
        sys.exit()
    ref=loadtxt('%s/target.txt' % dir.split('Data')[0])
    filename=genfile.split(dir)[1].split('.lh5')[0]
    if "Coarse" in filename:
        coarse=filename.split('Coarsed')[1].split('_')[0]
        rcut=filename.split('_r')[1].split('_')[0]
        modeldir='%s/msml%s_coarse_r%s_d%s/' % ( dir.split('/')[0], lag, rcut, coarse)
    else:
        modeldir='%s/msml%s' % (dir, lag)
    print "computing PMF from %s with lag %s frames" % (filename, lag)
    conf=Trajectory.load_from_pdb('%s/%s_noh.pdb' % (dir.split('Data')[0], system))
    lig_atoms=conf['XYZList'].shape[1]
    map=loadtxt('%s/Mapping.dat' % modeldir)
    pops=loadtxt('%s/Populations.dat' % modeldir)

    ligfile='%s.vmd_ligcoords.pickle' % genfile.split('.lh5')[0]
    if os.path.exists(ligfile):
        lighandle=open(ligfile, 'rb')
    else:
        print "need to run tcl scripts and pickler for coordinates"
        sys.exit()

    ligcoors=pickle.load(lighandle)
    lighandle.close()
    mapped_ligcoors, x_range, y_range, z_range, box_volume=get_minmax(ligcoors, map)
    correction=-0.6*log(box_volume/1600.0)
    # get prot-lig distances
    com_distances=loadtxt('%s.vmd_com.dat' % genfile.split('.lh5')[0], usecols=(1,))
    ref_com=com_distances[0]
    com_distances=com_distances[1:]
    if len(com_distances)==len(map):
        print "protein-ligand COM distances per state exist"
    # PMF
    frames=where(map!=-1)[0]
    mapped_states=map[frames]
    mapped_com=range(0, len(frames))
    mapped_com_distances=com_distances[frames]
    space=PMF3D.PMF3D(x_range, y_range, z_range)
    spacetrack=space.microstate_count_grid(mapped_ligcoors)
    new_pops={ key: pops[n] for (n, key) in enumerate(mapped_ligcoors.keys())}
    GD=space.new_manual_allcoor_grid(spacetrack, mapped_ligcoors, new_pops, type='pops')
    free=array([-0.6*log(i) for i in pops])
    subtract=min(free)
    free=array([k-subtract for k in free])
    GD=space.new_manual_allcoor_grid(spacetrack, mapped_ligcoors, new_pops, type='pops')
    GDfree=-0.6*log(GD)
    GDfree=PMF3D.convert(GDfree, max(free))
    GDfree=GDfree-min(GDfree.flatten())
    space.write_dx(GDfree, modeldir)
    frees=[]
    corrs=[]
    axis=[]
    volumes=[]
    cutoffs=arange(0, 40, 1)
    if volume=='all':
        for cutoff in cutoffs:
            bound_frames=where(mapped_com_distances < cutoff)[0]
            if len(bound_frames)==0:
                print "no bound states less than reference com distance %s" % cutoff
                continue
            gen_bound_frames=[]
            for i in bound_frames:
                location=where(map==i)[0]
                gen_bound_frames.append(location)
            new_points={ key: mapped_ligcoors[key] for key in bound_frames}
            new_pops={ key: pops[key] for key in bound_frames}
            GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
            boundspace=space.pmfvolume(GD)

            new_pops={ key: 1.0 for key in bound_frames}
            GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
            boundvolume=space.pmfvolume(GD)
            volumes.append(boundvolume)
            print "count bound volume ", boundvolume

            # unbound frames are above COM cutoff
            unbound_frames=array([int(x) for x in mapped_states if x not in bound_frames])
            new_points={ key: mapped_ligcoors[key] for key in unbound_frames}
            new_pops={ key: pops[key] for key in unbound_frames}
            GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
            unboundspace=space.pmfvolume(GD)

            # for counting states
            new_pops={ key: 1.0 for key in unbound_frames}
            GD=space.new_manual_allcoor_grid(spacetrack, new_points, new_pops, type='pops')
            unboundvolume=space.pmfvolume(GD)
            #print "unbound volume ", unboundvolume

            # free energy from ratio
            depth=-0.6*log(boundspace/unboundspace)
            frees.append(depth)
            axis.append(cutoff)
            print "corrected integrated dG ratio at cutoff %s is %s" % (cutoff, depth+correction)
        k=len(pops)
        savetxt('%s/%s_k%s_l%s_pocketdistance_free.dat' % (modeldir, filename, k, lag), [x+correction for x in frees])
        savetxt('%s/%s_k%s_l%s_pocketdistance_volume.dat' % (modeldir, filename, k, lag), volumes)
        savetxt('%s/%s_k%s_l%s_pocketdistance_axis.dat' % (modeldir, filename, k, lag), axis)
        pylab.figure()
        #pylab.plot(axis, corrs, color='green', label='correction')
        pylab.plot(axis, [ref]*len(frees), color='black', label='exp')
        pylab.plot(axis, [x+correction for x in frees], color='red', label='standard free')
        pylab.legend()
        pylab.savefig('%s/%s_k%s_l%s_pocketdistance_free.png' % (modeldir, filename, k, lag))
        #pylab.show()
    else:
        k=len(pops)
        standard_frees=loadtxt('%s/%s_k%s_l%s_pocketdistance_free.dat' % (modeldir, filename, k, lag))
        bound_volumes=loadtxt('%s/%s_k%s_l%s_pocketdistance_volume.dat' % (modeldir, filename, k, lag))
        axis=loadtxt('%s/%s_k%s_l%s_pocketdistance_axis.dat' % (modeldir, filename, k, lag))
        cutoff=float(volume)
        index=where(axis==cutoff)[0]
        print "standard free energy is %s" % (standard_frees[index])
        print "(population-weighted) bound volume is %s A^3" % bound_volumes[index]

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='system',
                      help='input system')
    parser.add_option('-g', '--genfile', dest='genfile',
                      help='input gens file lh5')
    parser.add_option('-l', '--lag', dest='lag',
                          help='lag time in steps')
    parser.add_option('-v', '--volume', dest='volume',
                          help='volume cutoff')
    parser.add_option('-p', '--nproc', dest='nproc',
                          help='num proc to run distance parallel')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(system=options.system, genfile=options.genfile, lag=options.lag, volume=options.volume, nproc=options.nproc)




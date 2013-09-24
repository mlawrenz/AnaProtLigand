from msmbuilder import io, Trajectory, Project
import sys
import os
import optparse
import numpy

def remap_ass(ref_ass, map_ass, unique, count):
    mapper=dict()
    for x in unique:
        if x==-1:
            pass
        else:
            locations=numpy.where(ref_ass==x)
            map_ass[locations]=count
            mapper[count]=x
            count+=1
    return map_ass, mapper, count

def make_map(index, all, new, map):
    return map

def main(coarse_val, orig_val, rcut):
    data=dict()
    data['coarse']=dict()
    data['orig']=dict()
    dirs=dict()
    dirs['coarse']='./d%s' % coarse_val
    dirs['orig']='./d%s' % orig_val
    proj=Project.load_from('ProjectInfo.yaml')
    types=['ass', 'rmsd', 'dist', 'gens']
    for key in ['coarse', 'orig']:
        for type in types:
            if 'ass' in type:
                ass=io.loadh('%s/Data/Assignments.h5' % dirs[key])
                data[key][type]=ass['arr_0']
            elif 'dist' in type:
                ass=io.loadh('%s/Data/Assignments.h5.distances' % dirs[key])
                data[key][type]=ass['arr_0']
            elif 'rmsd' in type:
                rmsd=numpy.loadtxt('%s/Gens.rmsd.dat' % dirs[key])
                data[key][type]=rmsd
            elif 'gens' in type:
                gens=Trajectory.load_from_lhdf('%s/Gens.lh5' % dirs[key])
                data[key][type]=gens
    unboundmap=dict()
    boundmap=dict()
    # build map dict for orig to coarse unbound states, bound will stay same
    unboundass=-1*numpy.ones(( data['orig']['ass'].shape[0], data['orig']['ass'].shape[1]), dtype=int)
    newass=-1*numpy.ones(( data['orig']['ass'].shape[0], data['orig']['ass'].shape[1]), dtype=int)
    newdist=-1*numpy.ones(( data['orig']['ass'].shape[0], data['orig']['ass'].shape[1]))
    for j in range(0, data['orig']['ass'].shape[0]):
        rmsd=numpy.loadtxt('Trajectories-metric/trj%s_lprmsd.dat' % j)
        frames=numpy.where(data['orig']['ass'][j]!=-1)[0]
        if len(rmsd)!=len(frames):
            print "trajectory mismatch"
            import pdb
            pdb.set_trace()
        for (n,i) in enumerate(data['orig']['ass'][j]):
            # if unbound
            if i != -1:
                #if data['orig']['rmsd'][i] > float(rcut):
                if rmsd[n] > float(rcut):
                    newstate=data['coarse']['ass'][j][n]
                    if data['coarse']['rmsd'][newstate] < float(rcut):
                        newass[j][n]=i
                        newdist[j][n]=data['orig']['dist'][j][n]
                    else:
                        unboundass[j][n]=newstate
                        newdist[j][n]=data['coarse']['dist'][j][n]
                else:
                    newass[j][n]=i
                    newdist[j][n]=data['orig']['dist'][j][n]
    count=0
    unique=sorted(set(newass.flatten()))
    newass, boundmap, count=remap_ass(newass, newass, unique, count)
    unique=sorted(set(unboundass.flatten()))
    newass, unboundmap, count=remap_ass(unboundass, newass, unique, count)
    import pdb
    pdb.set_trace()
    io.saveh('%s/Coarsed_r%s_d%s_Assignments_test.h5' % (dirs['orig'], rcut, coarse_val), newass)
    io.saveh('%s/Coarsed_r%s_d%s_Assignments_test.distances.h5' % (dirs['orig'], rcut, coarse_val), newdist)
    subdir='%s/Coarsed_r%s_gen/' % (dirs['orig'], rcut)
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    ohandle=open('%s/Coarsed%s_r%s_Gens_test.rmsd.dat' % (subdir, coarse_val, rcut), 'w')
    b=data['orig']['gens']['XYZList'].shape[1]
    c=data['orig']['gens']['XYZList'].shape[2]
    dicts=[boundmap, unboundmap]
    names=['bound', 'unbound']
    labels=['orig', 'coarse']
    total=len(boundmap.keys()) + len(unboundmap.keys())
    structure=proj.empty_traj()
    structure['XYZList']=numpy.zeros((total, b, c), dtype='float32')
    count=0
    for (name, label, mapdata) in zip( names, labels, dicts):
        print "writing coarse gen %s out of %s pdbs" % (count, len(mapdata.keys()))
        for i in sorted(mapdata.keys()):
            macro=mapdata[i]
            structure['XYZList'][count]=data[label]['gens']['XYZList'][macro]
            ohandle.write('%s\t%s\t%s\n' % (name, count, data[label]['rmsd'][macro]))
            print name, count
            count+=1
    otraj='%s/Coarsed%s_r%s_Gens_test.xtc' % (subdir, coarse_val, rcut)
    if os.path.exists(otraj):
        os.remove(otraj)
    structure.save_to_xtc('%s/Coarsed%s_r%s_Gens_test.xtc' % (subdir, coarse_val, rcut))
    #return data

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-r', '--rcut', dest='rcut',
                      help='unbound state distance cutoff')
    parser.add_option('-c', '--coarse', dest='coarse_val',
                      help='coarse clustering cutoff')
    parser.add_option('-o', '--orig', dest='orig_val',
                      help='original clustering cutoff')
    (options, args) = parser.parse_args()
    return (options, args)


if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(coarse_val=options.coarse_val, orig_val=options.orig_val, rcut=options.rcut)


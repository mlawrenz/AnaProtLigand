from msmbuilder import io, Trajectory, Project
import sys
import os
import optparse
import numpy

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
                data[key][type]=ass
            elif 'dist' in type:
                ass=io.loadh('%s/Data/Assignments.h5.distances' % dirs[key])
                data[key][type]=ass
            elif 'rmsd' in type:
                rmsd=numpy.loadtxt('%s/Gens.rmsd.dat' % dirs[key])
                data[key][type]=rmsd
            elif 'gens' in type:
                gens=Trajectory.load_from_lhdf('%s/Gens.lh5' % dirs[key])
                data[key][type]=gens
    unboundmap=dict()
    boundmap=dict()
    unboundstates=dict()
    unboundrmsd=dict()
    # build map dict for orig to coarse unbound states, bound will stay same
    for j in range(0, data['orig']['ass']['arr_0'].shape[0]):
        for (n,i) in enumerate(data['orig']['ass']['arr_0'][j]):
            if i!=-1:
                if data['orig']['rmsd'][i] > float(rcut):
                    state=data['coarse']['ass']['arr_0'][j][i]
                    if state!=-1:
                        if i not in unboundstates.keys():
                            unboundstates[i]=state
                            unboundstates[i]=data['coarse']['dist']['arr_0'][j][i]
                        elif i in unboundstates.keys():
                            if unboundstates[i]!=state:
                                prev=unboundrmsd[i]
                                if data['coarse']['dist']['arr_0'][j][i] < prev:
                                    unboundstates[i]=state
                                    unboundstates[i]=data['coarse']['dist']['arr_0'][j][i]
                                else:
                                    pass # keep previous assignment
                else:
                    if i not in boundmap.keys():
                        boundmap[i]=i
    import pdb
    pdb.set_trace()
    # deal with micro states that are sent to different coarse states
    map=dict()
    bound=0
    for i in sorted(boundmap.keys()):
        if i not in map.keys():
            map[i]=bound
            bound+=1
    unbound=max(map.values())+1
    unbound_start=max(map.values())+1
    unboundtrack=dict()
    for state in sorted(unboundmap.keys()):
        for i in unboundmap[state]:
            map[i]=unbound
            if i not in unboundtrack.keys():
                unboundtrack[i]=state
            elif i in unboundtrack.keys() and unboundtrack[i]!=state:
                print "state mismatch"
                import pdb
                pdb.set_trace()
        unbound+=1
    newgens=proj.empty_traj()
    subdir='%s/Coarsed_r%s_PDBs/' % (dirs['orig'], rcut)
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    ohandle=open('%s/Coarsed%s_r%s_Gens.rmsd.dat' % (subdir, coarse_val, rcut), 'w')
    b=data['orig']['gens']['XYZList'].shape[1]
    c=data['orig']['gens']['XYZList'].shape[2]
    ## add to newgens accoring to KEY values and their mapped states, will have multiple states per macrostate dict VALUE
    ## loop over values not keys, maybe write out PDBs of states?
    print "unbound start %s" % unbound_start
    for macro in sorted(set(map.values())):
        print "writing coarse gen %s out of %s pdbs" % (macro, len(set(map.values())))
        if macro >= unbound_start:
            import pdb
            pdb.set_trace()
            gens=proj.empty_traj()
            a=len(sorted(map.keys()))
            newrmsd=numpy.zeros((a))
            states=[i for i in map.keys() if map[i]==macro]
            structure=proj.empty_traj()
            structure['XYZList']=numpy.zeros((len(states), b, c), dtype=numpy.float32)
            for (n, state) in enumerate(states):
                unboundstate=unboundtrack[state]
                structure['XYZList'][n]=data['coarse']['gens']['XYZList'][unboundstate]
                structure.save_to_pdb('%s/Coarsed%s_r%s_Gen%s-%s.pdb' % (subdir, coarse_val, rcut, macro, n))
                ohandle.write('%s\t%s\t%s\n' % (macro, n, data['orig']['rmsd'][unboundstate]))
        else:
            states=[i for i in map.keys() if map[i]==macro]
            if len(states) > 1 :
                print "wrong"
                import pdb
                pdb.set_trace()
            state=states[0]
            structure=proj.empty_traj()
            structure['XYZList']=numpy.zeros((1, b, c), dtype=numpy.float32)
            structure['XYZList'][0]=data['orig']['gens']['XYZList'][state]
            structure.save_to_pdb('%s/Coarsed%s_r%s_Gen%s.pdb' % (subdir, coarse_val, rcut, macro))
            ohandle.write('%s\t%s\n' % (macro, data['orig']['rmsd'][state]))
    newass=numpy.zeros(( data['orig']['ass']['arr_0'].shape[0], data['orig']['ass']['arr_0'].shape[1]), dtype=int)
    for j in range(0, data['orig']['ass']['arr_0'].shape[0]):
        for (n,i) in enumerate(data['orig']['ass']['arr_0'][j]):
            if i!=-1:
                newass[j][n]=map[i]
    io.saveh('%s/Coarsed%s_r%s_Assignments.h5' % (subdir, rcut, coarse_val), newass)
    return data

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


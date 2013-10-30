import numpy, pylab
from msmbuilder import Trajectory, Project
import optparse

def main(dir, coarse, lag):
    data=dict()
    data['dist']=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.prot_lig_distance.dat' % (dir, coarse), usecols=(1,))
    rmsd=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.rmsd.dat' % (dir, coarse))
    data['rmsd']=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.rmsd.dat' % (dir, coarse))
    data['helixrmsd']=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.helixrmsd.dat' % (dir, coarse))
    com=numpy.loadtxt('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.vmd_com.dat' % (dir, coarse), usecols=(1,))
    com=[i/com[0] for i in com]
    data['com']=com[1:]
    pops=numpy.loadtxt('%s/msml%s_coarse_r10_d%s/Populations.dat' % (dir, lag, coarse))
    map=numpy.loadtxt('%s/msml%s_coarse_r10_d%s/Mapping.dat' % (dir, lag, coarse))
    
    data['com']=numpy.array(data['com'])
    data['rmsd']=numpy.array(data['rmsd'])
    data['rmsd']=numpy.array(data['rmsd'])
    #pylab.scatter(data['com'], data['rmsd'], c='k')  #pops, cmap=cm, s=50, alpha=0.7)
    #pylab.show()
    unbound_touch=numpy.where(data['dist'] > 4.0)[0]
    unbound_far=numpy.where(data['rmsd'] > 50.0)[0]
    pdb=Trajectory.load_from_pdb('sir2_bound_reference.pdb')
    proj=Project.load_from('ProjectInfo.yaml')
    count=0

    print unbound_far
    print "%s non touching unbound states" % count

    bound1=numpy.where((data['helixrmsd']<1.5))[0] # &(data['com']<5.0))[0]
    print "strict bound", bound1
    bound2=numpy.where((data['helixrmsd']<2.0))[0] # &(data['com']<5.0))[0]
    print "loose bound", bound2
    for i in bound1:
        if map[i]!=-1:
            count+=1
            t=Trajectory.load_from_xtc('%s/Coarsed_r10_gen/Coarsed%s_r10_Gens.xtc' % (dir, coarse), Conf=pdb)
            p=proj.empty_traj()
            p['XYZList']=numpy.zeros((1, t['XYZList'].shape[1], t['XYZList'].shape[2]), dtype=numpy.float32)
            p['XYZList'][0]=t['XYZList'][i]
            p.save_to_pdb('%s/bound_state%s.pdb' % (dir, int(map[i])))
    #print data['rmsd'][bound2]
    #print data['com'][bound2]
        
    cm=pylab.cm.get_cmap('RdYlBu_r') #blue will be negative components, red positive
    pylab.scatter(data['com'], data['rmsd'], c='k')  #pops, cmap=cm, s=50, alpha=0.7)
    pylab.scatter(data['com'][unbound_touch], data['rmsd'][unbound_touch], c='b')
    pylab.hold(True)
    #pylab.scatter(data['com'][bound1], data['rmsd'][bound1], c='purple')
    pylab.scatter(data['com'][bound2], data['rmsd'][bound2], c='r')
    pylab.ylabel('p53 RMSD wrt Protein')
    pylab.xlabel('p53-sir2 active site COM')
    #pylab.ylabel('p53 SELFRMSD')
   # 
    #pylab.figure()
    #pylab.scatter(data['com'], rmsd, c='k')  #pops, cmap=cm, s=50, alpha=0.7)
    #pylab.scatter(data['com'][unbound_touch], data['rmsd'][unbound_touch], c='b')
    #pylab.hold(True)
    #pylab.scatter(data['com'][bound2], data['rmsd'][bound2], c='r')
    #pylab.xlabel('p53-sir2 active site COM')
    pylab.show()
        
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    parser.add_option('-c', '--coarse', dest='coarse',
                      help='coarse grain cutoff')
    parser.add_option('-l', '--lag', dest='lag',
                      help='lag time')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(dir=options.dir, coarse=options.coarse, lag=options.lag)
                                                                                                        

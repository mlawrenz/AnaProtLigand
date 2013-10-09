#!/usr/bin/python
import pylab
from msmbuilder import io
import optparse
import numpy


def main(modeldir, gensfile, rcut=None):
    mapdata=dict()
    ass=io.loadh('%s/Assignments.Fixed.h5' % modeldir)

    data=dict()
    data['rmsd']=numpy.loadtxt('%s.rmsd.dat' % gensfile.split('.lh5')[0])
    com=numpy.loadtxt('%s.vmd_com.dat' % gensfile.split('.lh5')[0], usecols=(1,))
    data['com']=com[1:]
    pops=numpy.loadtxt('%s/Populations.dat' % modeldir)
    map=numpy.loadtxt('%s/Mapping.dat' % modeldir)

    mapdata['rmsd']=[]
    mapdata['com']=[]
    for x in range(0, len(data['rmsd'])):
        if map[x]!=-1:
            mapdata['com'].append(data['com'][x])
            mapdata['rmsd'].append(data['rmsd'][x])

    #RMSD cutoff 
    cutoffs=numpy.arange(1,30,0.5)
    bound_pops=[]
    for type in mapdata.keys(): 
        pylab.figure()
        ohandle=open('%s/%s_msm_frees.dat' % (modeldir, type), 'w')
        data=[]
        for cutoff in cutoffs:
            bound_pops=[]
            for (state, x) in enumerate(mapdata['rmsd']):
                if x < cutoff:
                    bound_pops.append(pops[state])
            ### calculate binding free energy from populations
            if len(bound_pops)==0:
                dG=100
            else:
                bound=numpy.sum(bound_pops)
                unbound=1-bound
                dG=-0.6*numpy.log(bound/unbound)
                #dG=-0.6*numpy.log(bound/(unbound**2))
            
            ### calculate standard state correction, in ansgtroms here
            boxvolume=244.80*(10**3)
            v0=1600
            corr=-0.6*numpy.log(boxvolume/v0)
            dG_corr=dG+corr 
            if cutoff==float(rcut):
                print cutoff, dG_corr
            data.append(dG_corr)
            ohandle.write('%s\t%s\n' % (cutoff, dG_corr))
        pylab.plot(cutoffs, data, label=type)
        pylab.legend()
        pylab.ylim(-8, (-1)*corr)
    pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--modeldir', dest='modeldir',
                      help='msm directory')
    parser.add_option('-g', '--gensfile', dest='gensfile',
                      help='gens files')
    parser.add_option('-c', '--cutoff', dest='cutoff',
                      help='cutoff')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(modeldir=options.modeldir, gensfile=options.gensfile, rcut=options.cutoff)


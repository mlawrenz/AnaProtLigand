from parallel_distance import parallel_distance
import os
import optparse
import pickle
import multiprocessing


def main(protfile, ligfile, nproc):
    dir=os.path.dirname(protfile)
    outname='%s.prot_lig_distance' % protfile.split('.vmd')[0]
    lighandle=open(ligfile, 'rb')
    ligcoors=pickle.load(lighandle)
    lighandle.close()
    prothandle=open(protfile, 'rb')
    protcoors=pickle.load(prothandle)
    prothandle.close()

    # get prot-lig distances
    statefile='%s.dat' % outname
    checkfile='%s.pickle.chkpt' % outname
    if os.path.exists(statefile) and os.path.getsize(statefile) > 0 :
        print "protein-ligand distances per state exist"
        states=loadtxt(statefile, usecols=(0,))
        state_distances=loadtxt(statefile, usecols=(1,))
    elif os.path.exists(checkfile) and os.path.getsize(checkfile) > 0 :
        ohandle=open(checkfile, 'rb')
        completed=pickle.load(ohandle)
        ohandle.close()
        print "restarting protein-ligand distances from checkpoint"
        completed=parallel_distance(nproc, ligcoors, protcoors, checkfile, completed)
        states=sorted(completed.keys())
        state_distances=[completed[i] for i in sorted(completed.keys())]
        statehandle=open(statefile, 'w')
        for state in sorted(completed.keys()):
            statehandle.write('%s\t%s\n' % (state, completed[state]))
    else:
        print "starting protein-ligand distances"
        completed=parallel_distance(nproc, ligcoors, protcoors, checkfile, completed=None)
        states=sorted(completed.keys())
        state_distances=[completed[i] for i in sorted(completed.keys())]
        statehandle=open(statefile, 'w')
        for state in sorted(completed.keys()):
            statehandle.write('%s\t%s\n' % (state, completed[state]))
    return completed

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-p', '--protfile', dest='protfile',
                      help='file with protein coors')
    parser.add_option('-l', '--ligfile', dest='ligfile',
                      help='file with ligand coors')
    parser.add_option('-n', '--nproc', dest='nproc',
                      help='n parallel processors')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(protfile=options.protfile, ligfile=options.ligfile, nproc=options.nproc)


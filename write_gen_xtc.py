from msmbuilder import Trajectory
import optparse, numpy

def main(dir):
    t=Trajectory.load_from_lhdf('%s/Gens.lh5' % dir)
    t.save_to_xtc('%s/Gens.xtc' % dir)

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    (options, args) = parser.parse_args()
    return (options, args)

#run the function main if namespace is main
if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(dir=options.dir)


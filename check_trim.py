import numpy
data=numpy.loadtxt('d6/Gens.rmsd.dat')
map=numpy.loadtxt('d6/msml1050_coarsed20_r10/Mapping.dat', dtype=int)
frames=numpy.where(map!=-1)[0]
tmp=numpy.where(data[frames]< 10)[0]
print "bind", data[frames][tmp]
tmp=numpy.where(data[frames]> 30)[0]
print "unbind", data[frames][tmp]
map=numpy.loadtxt('d6/msml1050_coarsed12_r10/Mapping.dat', dtype=int)
tmp=numpy.where(data[frames]< 10)[0]
print "bind", data[frames][tmp]
frames=numpy.where(map!=-1)[0]
tmp=numpy.where(data[frames]> 30)[0]
print "unbind", data[frames][tmp]

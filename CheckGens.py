import numpy, pylab
data=dict()
rmsd=numpy.loadtxt('d6/Coarsed_r10_gen/Coarsed20_r10_Gens.rmsd.dat', usecols=(2,))
data['rmsd']=numpy.loadtxt('d6/Coarsed_r10_gen/Coarsed20_r10_Gens.selfrmsd.dat')
com=numpy.loadtxt('d6/Coarsed_r10_gen/Coarsed20_r10_Gens.vmd_com.dat', usecols=(1,))
com=[i/com[0] for i in com]
data['com']=com[1:]
pops=numpy.loadtxt('d6/msml500_coarse_r10_d20/Populations.dat')
map=numpy.loadtxt('d6/msml500_coarse_r10_d20/Mapping.dat')

data['com']=numpy.array(data['com'])
data['rmsd']=numpy.array(data['rmsd'])
unbound1=numpy.where((data['rmsd']>7.0)&(data['com']>3.0))[0]
print "strict unbound", unbound1
unbound2=numpy.where((data['rmsd']>6.0)&(data['com']>2.0))[0]
print "loose unbound", unbound1
bound1=numpy.where((data['rmsd']<3.0)&(data['com']<1.1))[0]
print "strict bound", bound1
#print data['rmsd'][bound1]
#print data['com'][bound1]
bound2=numpy.where((data['rmsd']<4.0)&(data['com']<1.1))[0]
print "loose bound", bound2
#print data['rmsd'][bound2]
#print data['com'][bound2]

cm=pylab.cm.get_cmap('RdYlBu_r') #blue will be negative components, red positive
pylab.scatter(data['com'], data['rmsd'], c='k')  #pops, cmap=cm, s=50, alpha=0.7)
pylab.scatter(data['com'][unbound2], data['rmsd'][unbound2], c='b')
pylab.hold(True)
#pylab.scatter(data['com'][bound1], data['rmsd'][bound1], c='purple')
pylab.scatter(data['com'][bound2], data['rmsd'][bound2], c='r')
pylab.xlabel('p53-s100b active site COM')
pylab.ylabel('p53 RMSD')

pylab.figure()
pylab.scatter(data['com'], rmsd, c='k')  #pops, cmap=cm, s=50, alpha=0.7)
pylab.scatter(data['com'][unbound2], rmsd[unbound2], c='b')
pylab.hold(True)
pylab.scatter(data['com'][bound2], rmsd[bound2], c='r')
pylab.xlabel('p53-s100b active site COM')
pylab.ylabel('p53 RMSD wrt Protein')
pylab.show()
        

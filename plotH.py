import numpy, pylab
import scipy.ndimage as ndimage
rmsd=numpy.loadtxt('all_rmsd')
dist=numpy.loadtxt('all_dist')
#bbrmsd=numpy.loadtxt('trj-rmsd/all_bbrmsd')

fig=pylab.figure()
halfdist=dist[0:int(len(dist)/2)]
halfrmsd=rmsd[0:int(len(rmsd)/2)]
xbins=numpy.arange(5.5, 43, 2)
ybins=numpy.arange(0.5, 56, 2)
(H, x, y)=numpy.histogram2d(halfdist, halfrmsd,  bins=[numpy.arange(5.5, 43, 2), numpy.arange(0.5, 56, 2)], normed=True)
Hhalf=-0.6*numpy.log(H)
Hhalf=Hhalf-min(Hhalf.flatten())
#pylab.pcolor(Hnew.transpose())
#H2 = ndimage.gaussian_filter(Hnew, sigma=0.2, order=0)
#pylab.imshow(H2[:][::-1], interpolation='bilinear')#, interpolation='nearest')
## FULL DATA
(H, x, y)=numpy.histogram2d(dist, rmsd, bins=[numpy.arange(5.5, 43, 2), numpy.arange(0.5, 56, 2)], normed=True)
Hnew=-0.6*numpy.log(H)
Hnew=Hnew-min(Hnew.flatten())
index=numpy.where(numpy.array(sorted(Hhalf.flatten())[::-1])!=numpy.inf)[0][0]
Hhalfmax=sorted(Hhalf.flatten())[::-1][index]
index=numpy.where(numpy.array(sorted(Hnew.flatten())[::-1])!=numpy.inf)[0][0]
Hmax=sorted(Hnew.flatten())[::-1][index]
Hmin=min(Hnew.flatten())
if Hhalfmax>Hmax:
    Hmax=Hhalfmax
if min(Hhalf.flatten())<Hmin:
    Hmin=min(Hhalf.flatten())
print Hmax
print Hmin
where=numpy.isinf(Hhalf)
Hhalf[where]=Hmax
Hhalf=Hhalf.transpose()
where=numpy.isinf(Hnew)
Hnew[where]=Hmax
Hnew=Hnew.transpose()
ax=fig.add_subplot(1,2,1)
ax.imshow(Hhalf[:][::-1],vmin=Hmin, vmax=Hmax,  interpolation='nearest')
loc,labels=pylab.yticks(numpy.arange(0, len(ybins),4), ybins[::-1][::4])
loc,labels=pylab.xticks(numpy.arange(0, len(xbins),4), xbins[::4])
pylab.xlabel('(Half Data) Ligand N - F99 C$\\alpha$ Distance ($\AA$)')
pylab.ylabel('Ligand RMSD from Crystal ($\AA$)')
ax=fig.add_subplot(1,2,2)
#pylab.show()
other=ax.imshow(Hnew[:][::-1], vmin=Hmin, vmax=Hmax, interpolation='nearest')
loc,labels=pylab.yticks(numpy.arange(0, len(ybins),4), ybins[::-1][::4])
loc,labels=pylab.xticks(numpy.arange(0, len(xbins),4), xbins[::4])
pylab.xlabel('Ligand OH - F99 C$\\alpha$ Distance ($\AA$)')
#pylab.colorbar(other)
pylab.savefig('LG2_landscape_compare1.png', dpi=300)
#pylab.savefig('LG2_landscape_compare1_colorbar.png')
pylab.show()


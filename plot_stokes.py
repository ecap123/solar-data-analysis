# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import idlsave

 
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

figname = "180624_gregor_m1_stokesprofiles.eps"

s = idlsave.read('4227_m1_ff_ct_shift.dat') 
F=0.00412634
wl0=4223.95       
wlend=wl0+1239*F                      #wl corrispondente al pixel 1240   
wlarray=np.arange(wl0, wlend, F)

#averaged slit in some range

Msiy=(np.mean(s.si[110:130,0:1240],axis=0))  
Msiymax=np.max(Msiy)        #normalizzata per il continuo
emptyarray= np.empty(1240)   #creo un vettore vuoto di 1240 spazi
emptyarray.fill(Msiymax)     # E lo riempio con il valore massimo dell'intensita'

Msiynorm=Msiy/emptyarray     #normalizzazione dell'intensita': I/Imax

Msqy=(np.mean(s.sq[110:130,0:1240],axis=0))*100
Msuy=(np.mean(s.su[110:130,0:1240],axis=0))*100   
Msvy=(np.mean(s.sv[110:130,0:1240],axis=0))*100

cols =  "magenta"

#
strstr = r"$Bcore = -14 $G"
strstr1 = r"$Bwings = -13 $G"


#plot all stokes profiles

fig = plt.figure(figsize=(36,29.5))
plt.subplots_adjust(left=0.12,right=0.98,bottom=0.12,top=0.90,wspace=0.32,hspace=0.25)
fig.suptitle('$\mu = 0.09$',fontsize=70)

#two panels
#I
plt.subplot(221)
#plt.yscale('log')
plt.yscale('linear')
axi = plt.gca()
axi.get_xaxis().get_major_formatter().set_useOffset(False)
axi.tick_params(axis='x',which='major',pad=26,labelsize=52,length=11,width=3.75)
axi.tick_params(axis='y',which='major',pad=25,labelsize=52,length=11,width=3.75)
#t = axi.yaxis.get_offset_text()
#t.set_size(44)

for axis in ['top','bottom','left','right']:
 axi.spines[axis].set_linewidth(4.5)

#V
plt.subplot(222)
plt.yscale('linear')

axv = plt.gca()
axv.get_xaxis().get_major_formatter().set_useOffset(False)
axv.tick_params(axis='x',which='major',pad=26,labelsize=52,length=11,width=3.75)
axv.tick_params(axis='y',which='major',pad=25,labelsize=52,length=11,width=3.75)
#t = axv.yaxis.get_offset_text()
#t.set_size(24)

for axis in ['top','bottom','left','right']:
 axv.spines[axis].set_linewidth(4.5)

#Q
plt.subplot(223)
plt.yscale('linear')

axq = plt.gca()
axq.get_xaxis().get_major_formatter().set_useOffset(False)
axq.tick_params(axis='x',which='major',pad=26,labelsize=52,length=11,width=3.75)
axq.tick_params(axis='y',which='major',pad=25,labelsize=52,length=11,width=3.75)
#t = axq.yaxis.get_offset_text()
#t.set_size(44)

for axis in ['top','bottom','left','right']:
 axq.spines[axis].set_linewidth(4.5)

#U
plt.subplot(224)
plt.yscale('linear')

axu = plt.gca()
axu.get_xaxis().get_major_formatter().set_useOffset(False)
axu.tick_params(axis='x',which='major',pad=26,labelsize=52,length=11,width=3.75)
axu.tick_params(axis='y',which='major',pad=25,labelsize=52,length=11,width=3.75)
#t = axu.yaxis.get_offset_text()
#t.set_size(44)

for axis in ['top','bottom','left','right']:
 axu.spines[axis].set_linewidth(4.5)

#
axi.plot(wlarray,Msiynorm,lw=3.2,color=cols)
axq.plot(wlarray,Msqy,lw=3.2,color=cols)
axu.plot(wlarray,Msuy,lw=3.2,color=cols)
axv.plot(wlarray,Msvy,lw=3.2,color=cols)

#l'unità di misura nel nostro caso sono counts
axi.set_ylabel(r"$I/I_c$",fontsize=52) 
axi.text(0.035,0.18,strstr,horizontalalignment="left",verticalalignment='center',\
 transform=axi.transAxes,fontsize=50)
axi.text(0.025,0.1,strstr1,horizontalalignment="left",verticalalignment='center',\
 transform=axi.transAxes,fontsize=50)
axv.set_ylabel(r"$V/I$ (\%)",fontsize=52)
axq.set_ylabel(r"$Q/I$ (\%)",fontsize=52)
axq.set_xlabel(r"Wavelength (\AA)",fontsize=50)
axu.set_ylabel(r"$U/I$ (\%)",fontsize=52)
axu.set_xlabel(r"Wavelength (\AA)",fontsize=50)


axi.set_xlim(wlarray[180],wlarray[1150])
#axi.set_xlim(4224.7,4228.7)
#axi.set_ylim(0.8e3,1.e4)
axi.xaxis.set_major_locator(plt.MaxNLocator(5))
axv.set_xlim(wlarray[180],wlarray[1150])
axv.xaxis.set_major_locator(plt.MaxNLocator(5))
axq.set_xlim(wlarray[180],wlarray[1150])
axq.set_ylim(-0.5,0.85)
axq.xaxis.set_major_locator(plt.MaxNLocator(5))
axu.set_xlim(wlarray[180],wlarray[1150])
#axu.set_ylim(-5.e-6,0.004)
axu.xaxis.set_major_locator(plt.MaxNLocator(5))
fig.savefig(figname)
#fig.savefig(fname,format='eps',dpi=1000)
#plt.close(fig)
#plt.show()

#!/usr/bin/python

#plot dei range di ali selezionate:


import numpy as np
import os, fnmatch
import matplotlib.pyplot as plt
import astropy
import idlsave
import scipy.io as spio
import matplotlib.ticker as ticker
from scipy import stats
from matplotlib import rc
from PIL import Image

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\boldmath']

#s = idlsave.read('4227_m4_ff_shift.dat')  
s = idlsave.read('4227_m1_ff_ct_shift.dat')  
s.si.shape  #(140, 1240)= (num immagini, y, x)
s.sq.shape  #(140, 1240)
s.su.shape  #(140, 1240)
s.sv.shape  #(140, 1240)

F=0.00412634
wl0=4223.95      
wlend=wl0+1239*F                      #wl corrispondente al pixel 1240   
wlarray=np.arange(wl0, wlend, F)

cols =  "black"

#plot image all stokes

plt.rc('axes.formatter', useoffset=False)
fig = plt.figure(figsize=(36,22.25))#golden number=1.618 --> 36/1.618=22.25

plt.subplot(4,1,1)
plt.subplots_adjust(left = 0.11, bottom=0.08, right = 0.92, top = 0.94, wspace=None, hspace = 0.4)

plt.gray()
#gregor:1pixel=0.33", 140 pixel=46"
plt.imshow(s.si, aspect='auto', extent=[wl0,wlend,46.2,0.3])
plt.ylim(plt.ylim()[::-1])
plt.tick_params(axis='x',which='major',pad=7,labelsize=30,length=6,width=3)
plt.tick_params(axis='y',which='major',pad=7,labelsize=30,length=6,width=3)
plt.ylabel(r'\textbf{$I$ (")}',fontsize=30)
plt.xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=30)
#plt.title(r'\textbf{$6"$}',fontsize=35)  #opposite solar limb forse 6"


plt.subplot(4,1,2)
plt.gray()
plt.imshow(s.sq, aspect='auto', extent=[wl0,wlend,46.2,0.3], vmin=-0.012, vmax=0.023)
plt.ylim(plt.ylim()[::-1])
plt.ylabel(r'\textbf{$Q/I$ (")}',fontsize=30)
plt.xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=30)
plt.tick_params(axis='x',which='major',pad=7,labelsize=30,length=6,width=3)
plt.tick_params(axis='y',which='major',pad=7,labelsize=30,length=6,width=3)


plt.subplot(4,1,3)
plt.gray()
plt.imshow(s.su, aspect='auto', extent=[wl0,wlend,46.2,0.3], vmin=-0.012, vmax=0.023)
plt.ylim(plt.ylim()[::-1])
plt.ylabel(r'\textbf{$U/I$ (")}',fontsize=30)
plt.xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=30)
plt.tick_params(axis='x',which='major',pad=7,labelsize=30,length=6,width=3)
plt.tick_params(axis='y',which='major',pad=7,labelsize=30,length=6,width=3)

plt.subplot(4,1,4)
plt.gray()
plt.imshow(s.sv, aspect='auto', extent=[wl0,wlend,46.2,0.3], vmin=-0.012, vmax=0.023)
plt.ylim(plt.ylim()[::-1])
plt.ylabel(r'\textbf{$V/I$ (")}',fontsize=30)
plt.xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=30)
plt.tick_params(axis='x',which='major',pad=7,labelsize=30,length=6,width=3)
plt.tick_params(axis='y',which='major',pad=7,labelsize=30,length=6,width=3)


plt.savefig('4227_m1_stokimage.eps')
#plt.show()


##########################################


#all averaged slit (stokes profile)
Msiy=np.mean(s.si, axis=0)
Msiymax=np.max(Msiy)        #normalizzata per il continuo
emptyarray= np.empty(1240)   #creo un vettore vuoto di 1240 spazi
emptyarray.fill(Msiymax)     # E lo riempio con il valore massimo dell'intensita'

Msiynorm=Msiy/emptyarray     #normalizzazione dell'intensita': I/Imax

Msqy=np.mean(s.sq, axis=0)*100
Msuy=np.mean(s.su, axis=0)*100    #Msqy.shape (1240,) mediato su tutta la slit
Msvy=np.mean(s.sv, axis=0)*100


#
#fig = plt.figure(figsize=(36,29.5))
fig = plt.figure(figsize=(36,22.25))#golden number=1.618 --> 36/1.618=22.25
plt.subplots_adjust(left=0.095,right=0.96,bottom=0.21,top=0.91,wspace=0.4,hspace=0.3)
#The two panels
plt.subplot(221)
plt.yscale('linear')
axi = plt.gca()
axi.get_xaxis().get_major_formatter().set_useOffset(False)
axi.tick_params(axis='x',which='major',pad=26,labelsize=40,length=11,width=3.75)
axi.tick_params(axis='y',which='major',pad=25,labelsize=40,length=11,width=3.75)
t = axi.yaxis.get_offset_text()
t.set_size(50)
for axis in ['top','bottom','left','right']:
 axi.spines[axis].set_linewidth(3.5)
#Stokes V: Check, but likely best not to use it
plt.subplot(222)
plt.yscale('linear')
axv = plt.gca()
axv.get_xaxis().get_major_formatter().set_useOffset(False)
axv.tick_params(axis='x',which='major',pad=26,labelsize=40,length=11,width=3.75)
axv.tick_params(axis='y',which='major',pad=25,labelsize=40,length=11,width=3.75)
axv.ticklabel_format(axis='y',style='sci',scilimits=(-3,3))
t = axv.yaxis.get_offset_text()
t.set_size(50)
for axis in ['top','bottom','left','right']:
 axv.spines[axis].set_linewidth(3.5)
#Stokes Q
plt.subplot(223)
plt.yscale('linear')
axq = plt.gca()
axq.get_xaxis().get_major_formatter().set_useOffset(False)
axq.tick_params(axis='x',which='major',pad=26,labelsize=40,length=11,width=3.75)
axq.tick_params(axis='y',which='major',pad=25,labelsize=40,length=11,width=3.75)
axq.ticklabel_format(axis='y',style='sci',scilimits=(-3,3))
t = axq.yaxis.get_offset_text()
t.set_size(50)
for axis in ['top','bottom','left','right']:
 axq.spines[axis].set_linewidth(3.5)
#Stokes U
plt.subplot(224)
plt.yscale('linear')
axu = plt.gca()
axu.get_xaxis().get_major_formatter().set_useOffset(False)
axu.tick_params(axis='x',which='major',pad=26,labelsize=40,length=11,width=3.75)
axu.tick_params(axis='y',which='major',pad=25,labelsize=40,length=11,width=3.75)
axu.ticklabel_format(axis='y',style='sci',scilimits=(-3,3))
t = axu.yaxis.get_offset_text()
t.set_size(50)
for axis in ['top','bottom','left','right']:
 axu.spines[axis].set_linewidth(3.5)

#
axi.plot(wlarray,Msiynorm,lw=2.2,color=cols)

axq.plot(wlarray,Msqy,lw=2.2,color=cols)

axu.plot(wlarray,Msuy,lw=2.2,color=cols)

axv.plot(wlarray,Msvy,lw=2.2,color=cols)

#axi.set_title('6"')
axi.set_ylabel(r"\textbf{$<I/I_c>$} ",fontsize=40)
axi.set_xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=38)
axi.set_xlim(wl0,wlend)
#axi.set_ylim(2e-7,2.e-5)
axq.set_ylabel(r"\textbf{$<Q/I>$ (\%)}",fontsize=40)
axq.set_xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=38)
axq.set_xlim(wl0,wlend)
axq.set_ylim(-0.74,2.5)
axu.set_ylabel(r"\textbf{$<U/I>$ (\%)}",fontsize=40)
axu.set_xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=38)
axu.set_xlim(wl0,wlend)
axu.set_ylim(-0.74,2.5)
axv.set_ylabel(r"\textbf{$<V/I>$ (\%)}",fontsize=40)
axv.set_xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=38)
axv.set_xlim(wl0,wlend)
axv.set_ylim(-0.74,2.5)

plt.savefig('4227_m1_stokprofile.eps')

#plt.show()



#!/usr/bin/python
#frazione di polarizzazione lineare puntuale lungo la slit e frazione di polarizzazione lineare media lungo la slit

import numpy as np
import os, fnmatch
import matplotlib.pyplot as plt
import astropy
import idlsave
import scipy.io as spio
import matplotlib.ticker as ticker
from scipy import stats
import pylab
from matplotlib import rc
from PIL import Image

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\boldmath']


s = idlsave.read('4227_m1_ff_ct_shift.dat')  
s.si.shape  #(140, 1240)= (num immagini, y, x)
s.sq.shape  #(140, 1240)
s.su.shape  #(140, 1240)
s.sv.shape  #(140, 1240)

F=0.00412634
wl0=4223.95      
wlend=wl0+1239*F                      #wl corrispondente al pixel 1240   
wlarray=np.arange(wl0, wlend, F)


#Blue wing I : 
a1=508
b1=581
#Blue wing II :
a11=446
b11=480
#Blue wing III : 
a111=391
b111=425
#Red wing  : 
a2=707
b2=806
#Core : 
ac=664
bc=694



###########################################################
#frazione di pol. mediata lungo tutta la slit

plt.rc('axes.formatter', useoffset=False)
fig = plt.figure(figsize=(36,22.25))#golden number=1.618 --> 36/1.618=22.25


Msqy=np.mean(s.sq, axis=0)
Msuy=np.mean(s.su, axis=0)    #Msqy.shape (1240,)

#Blue wing I : 

polav=np.sqrt(np.square(Msqy[a1])+np.square(Msuy[a1]))
for x in range(a1+1,b1+1):  
	polav=np.append(polav,np.sqrt(np.square(Msqy[x])+np.square(Msuy[x])))

plt.subplot(311)
ax = plt.gca()
plt.subplots_adjust(left  = 0.13, bottom=None, right = 0.92, top = 0.94, wspace=None, hspace = 0.6)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)

ax.scatter(wlarray[a1:b1+1],polav*100, label='Blue wing I',s=200) 
#parafit=np.polyfit(wlarray[a1:b1+1],polav*100,1) #1=grado del polinomio (retta), da iparametri fit lineare (m,q).
#print(parafit[0])  #coeff angolare retta
#linfit=np.poly1d(parafit) #costruisce il fit, noti i parametri (m,q)
#ax.plot(wlarray[a1:b1+1],linfit(wlarray[a1:b1+1]),color='black')
ax.set_ylabel(r'\textbf{$<P_L>$\ (\%)}',fontsize=30)
ax.set_xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=30)
ax.set_title('4"')	
ax.legend(fontsize=47.5)	

#Red wing  : 

polav=np.sqrt(np.square(Msqy[a2])+np.square(Msuy[a2]))
for x in range(a2+1,b2+1):  
	polav=np.append(polav,np.sqrt(np.square(Msqy[x])+np.square(Msuy[x])))

plt.subplot(312)
ax = plt.gca()
plt.subplots_adjust(left  = 0.13, bottom=None, right = 0.92, top = 0.94, wspace=None, hspace = 0.6)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.scatter(wlarray[a2:b2+1],polav*100, label='Red wing ',s=200) 
#parafit=np.polyfit(wlarray[a2:b2+1],polav*100,1) 
#print(parafit[0])  #coeff angolare retta
#linfit=np.poly1d(parafit) #costruisce il fit, noti i parametri (m,q)
#ax.plot(wlarray[a2:b2+1],linfit(wlarray[a2:b2+1]),color='black')
ax.set_ylabel(r'\textbf{$<P_L>$\ (\%)}',fontsize=30)
ax.set_xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=30)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.legend(fontsize=47.5)	


#Core : 

polav=np.sqrt(np.square(Msqy[ac])+np.square(Msuy[ac]))
for x in range(ac+1,bc+1):  
	polav=np.append(polav,np.sqrt(np.square(Msqy[x])+np.square(Msuy[x])))

plt.subplot(313)
ax = plt.gca()
plt.subplots_adjust(left  = 0.13, bottom=None, right = 0.92, top = 0.94, wspace=None, hspace = 0.6)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.scatter(wlarray[ac:bc+1],polav*100, label='Core',s=200) 
#parafit=np.polyfit(wlarray[ac:bc+1],polav*100,1) 
#print(parafit[0])  #coeff angolare retta
#linfit=np.poly1d(parafit) #costruisce il fit, noti i parametri (m,q)
#ax.plot(wlarray[ac:bc+1],linfit(wlarray[ac:bc+1]),color='black')
ax.set_ylabel(r'\textbf{$<P_L>$\ (\%)}',fontsize=30)
ax.set_xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=30)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.legend(fontsize=47.5)

plt.savefig('4227_m1_pol_all_av_slit.png')
#plt.show() 

###########################################################
#frazione di pol. mediata lungo tutta la slit per le Blue wings
plt.rc('axes.formatter', useoffset=False)
fig = plt.figure(figsize=(36,22.25))#golden number=1.618 --> 36/1.618=22.25


Msqy=np.mean(s.sq, axis=0)
Msuy=np.mean(s.su, axis=0)    #Msqy.shape (1240,)

#Blue wing I : 
polav=np.sqrt(np.square(Msqy[a1])+np.square(Msuy[a1]))
for x in range(a1+1,b1+1):  
	polav=np.append(polav,np.sqrt(np.square(Msqy[x])+np.square(Msuy[x])))


plt.subplot(311)
ax = plt.gca()
plt.subplots_adjust(left  = 0.13, bottom=None, right = 0.92, top = 0.94, wspace=None, hspace = 0.6)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.scatter(wlarray[a1:b1+1],polav*100, label='Blue wing I',s=200) 
#parafit=np.polyfit(wlarray[a1:b1+1],polav*100,1) #1=grado del polinomio (retta), da iparametri fit lineare (m,q).
#print(parafit[0])  #coeff angolare retta
#linfit=np.poly1d(parafit) #costruisce il fit, noti i parametri (m,q)
#ax.plot(wlarray[a1:b1+1],linfit(wlarray[a1:b1+1]),color='black')
ax.set_ylabel(r'\textbf{$<P_L>$\ (\%)}',fontsize=30)
ax.set_xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=30)
ax.set_title('4"')	
ax.legend(fontsize=47.5)	


#Blue wing II

polav=np.sqrt(np.square(Msqy[a11])+np.square(Msuy[a11]))
for x in range(a11+1,b11+1):  
	polav=np.append(polav,np.sqrt(np.square(Msqy[x])+np.square(Msuy[x])))

plt.subplot(312)
ax = plt.gca()
plt.subplots_adjust(left  = 0.13, bottom=None, right = 0.92, top = 0.94, wspace=None, hspace = 0.6)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.scatter(wlarray[a11:b11+1],polav*100, label='Blue wing II',s=200) 
#parafit=np.polyfit(wlarray[a1:b1+1],polav*100,1) #1=grado del polinomio (retta), da iparametri fit lineare (m,q).
#print(parafit[0])  #coeff angolare retta
#linfit=np.poly1d(parafit) #costruisce il fit, noti i parametri (m,q)
#ax.plot(wlarray[a1:b1+1],linfit(wlarray[a1:b1+1]),color='black')
ax.set_ylabel(r'\textbf{$<P_L>$\ (\%)}',fontsize=30)
ax.set_xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=30)
ax.legend(fontsize=47.5)	



#Blue wing III

polav=np.sqrt(np.square(Msqy[a111])+np.square(Msuy[a111]))
for x in range(a111+1,b111+1):  
	polav=np.append(polav,np.sqrt(np.square(Msqy[x])+np.square(Msuy[x])))

plt.subplot(313)
ax = plt.gca()
plt.subplots_adjust(left  = 0.13, bottom=None, right = 0.92, top = 0.94, wspace=None, hspace = 0.6)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.scatter(wlarray[a111:b111+1],polav*100, label='Blue wing III',s=200) 
#parafit=np.polyfit(wlarray[a1:b1+1],polav*100,1) #1=grado del polinomio (retta), da iparametri fit lineare (m,q).
#print(parafit[0])  #coeff angolare retta
#linfit=np.poly1d(parafit) #costruisce il fit, noti i parametri (m,q)
#ax.plot(wlarray[a1:b1+1],linfit(wlarray[a1:b1+1]),color='black')
ax.set_ylabel(r'\textbf{$<P_L>$\ (\%)}',fontsize=30)
ax.set_xlabel(r"\textbf{$Wavelength$\ (\AA)}",fontsize=30)
ax.legend(fontsize=47.5)	


plt.savefig('4227_m1_pol_all_av_slit_tris.png')
#plt.show() 

###########################################################
###########################################################
#frazone di pol. calcolata mediando lungo le wl e lasciando correre lungo la slit


plt.rc('axes.formatter', useoffset=False)
fig = plt.figure(figsize=(36,22.25))#golden number=1.618 --> 36/1.618=22.25

#Blue wing I 
sqq1=s.sq[0:140, a1:b1+1]   
Msqx1=np.mean(sqq1, axis=1) # medio lungo le x (colonne)
suu1=s.su[0:140, a1:b1+1]
Msux1=np.mean(suu1, axis=1) 


polavb=np.sqrt(np.square(Msqx1[0])+np.square(Msux1[0]))
for y in range(0,139):  
	polavb=np.append(polavb,np.sqrt(np.square(Msqx1[y])+np.square(Msux1[y])))
print polavb.shape
#gregor:1pixel=0.33", 140 pixel=46.2"  :

arcsecond=np.arange(0.3,46.2,0.33)

ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.plot(arcsecond,polavb*100, label='Blue wing I',color='blue')

#ax.set_ylabel(r'\textbf{$P_L$\ (\%)}',fontsize=30)
#ax.set_xlabel(r'$Slit\ position$\ (")',fontsize=30)
ax.legend(loc='upper center',fontsize=47.5) 
#ax.set_title('4"')


#Red wing 
sqq2=s.sq[0:140, a2:b2+1]   
Msqx2=np.mean(sqq2, axis=1) 
suu2=s.su[0:140, a2:b2+1]
Msux2=np.mean(suu2, axis=1)

polavr=np.sqrt(np.square(Msqx2[0])+np.square(Msux2[0]))
for y in range(0,139):  
	polavr=np.append(polavr,np.sqrt(np.square(Msqx2[y])+np.square(Msux2[y])))


ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.plot(arcsecond,polavr*100, label='Red wing',color='red')

#ax.set_ylabel(r'\textbf{$P_L$\ (\%)}',fontsize=30)
#ax.set_xlabel(r'$Slit\ position$\ (")',fontsize=30)
ax.legend(loc='upper center',fontsize=47.5) 
#ax.set_title('4"')


#Core : 
sqqc=s.sq[0:140, ac:bc+1]   
Msqxc=np.mean(sqqc, axis=1) 
suuc=s.su[0:140, ac:bc+1]
Msuxc=np.mean(suuc, axis=1)

polavc=np.sqrt(np.square(Msqxc[0])+np.square(Msuxc[0]))
for y in range(0,139):  
	polavc=np.append(polavc,np.sqrt(np.square(Msqxc[y])+np.square(Msuxc[y])))


ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.plot(arcsecond,polavc*100, label='Core', color='black')

ax.set_ylabel(r'\textbf{$P_L$\ (\%)}',fontsize=30)
ax.set_xlabel(r'$Slit\ position$\ (")',fontsize=30)
ax.legend(loc='upper right',fontsize=47.5) 
ax.set_xlim(arcsecond[0],arcsecond[139])
#ax.set_title('4"')


plt.savefig('4227_m1_pol_wl_averaged.png')

###########################################################
###########################################################
#frazone di pol. calcolata mediando lungo le wl e lasciando correre lungo la slit per le Blue wings

plt.rc('axes.formatter', useoffset=False)
fig = plt.figure(figsize=(36,22.25))#golden number=1.618 --> 36/1.618=22.25

#Blue wing I 
sqq1=s.sq[0:140, a1:b1+1]   
Msqx1=np.mean(sqq1, axis=1) # medio lungo le x (colonne)
suu1=s.su[0:140, a1:b1+1]
Msux1=np.mean(suu1, axis=1) 

polavb1=np.sqrt(np.square(Msqx1[0])+np.square(Msux1[0]))
for y in range(0,139):  
	polavb1=np.append(polavb1,np.sqrt(np.square(Msqx1[y])+np.square(Msux1[y])))

arcsecond=np.arange(0.3,46.2,0.33)

ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.plot(arcsecond,polavb1*100, label='Blue wing I',color='blue')

#ax.set_ylabel(r'\textbf{$P_L$\ (\%)}',fontsize=30)
#ax.set_xlabel(r'$Slit\ position$\ (")',fontsize=30)
#ax.set_title('4"')



#Blue wing II
sqq11=s.sq[0:140, a11:b11+1]   
Msqx11=np.mean(sqq2, axis=1) 
suu11=s.su[0:140, a11:b11+1]
Msux11=np.mean(suu11, axis=1)

polavb2=np.sqrt(np.square(Msqx11[0])+np.square(Msux11[0]))
for y in range(0,139):  
	polavb2=np.append(polavb2,np.sqrt(np.square(Msqx11[y])+np.square(Msux11[y])))


ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.plot(arcsecond,polavb2*100, label='Blue wing II',color='gray')

#ax.set_ylabel(r'\textbf{$P_L$\ (\%)}',fontsize=30)
#ax.set_xlabel(r'$Slit\ position$\ (")',fontsize=30)
#ax.set_title('4"')




#Blue wing III
sqq111=s.sq[0:140, a111:b111+1]   
Msqx111=np.mean(sqq111, axis=1) 
suu111=s.su[0:140, a111:b111+1]
Msux111=np.mean(suu111, axis=1)

polavb3=np.sqrt(np.square(Msqx111[0])+np.square(Msux111[0]))
for y in range(0,139):  
	polavb3=np.append(polavb3,np.sqrt(np.square(Msqx111[y])+np.square(Msux111[y])))

ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.tick_params(axis='x',which='major',pad=26,labelsize=45,length=11,width=3.75)
ax.tick_params(axis='y',which='major',pad=25,labelsize=45,length=11,width=3.75)
t = ax.yaxis.get_offset_text()
t.set_size(45)
for axis in ['top','bottom','left','right']:
 ax.spines[axis].set_linewidth(5.5)


ax.plot(arcsecond,polavb3*100, label='Blue wing III', color='green')

ax.set_ylabel(r'\textbf{$P_L$\ (\%)}',fontsize=30)
ax.set_xlabel(r'$Slit\ position$\ (")',fontsize=30)
ax.legend(loc='upper right',fontsize=47.5) 
ax.set_xlim(arcsecond[0],arcsecond[139])
#ax.set_title('4"')


plt.savefig('4227_m1_pol_wl_averaged_tris.png')
plt.show()

 
 



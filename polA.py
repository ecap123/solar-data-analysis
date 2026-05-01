#!/usr/bin/python

#angolo di polarizzazione mediando q e u lungo tutta la slit e lasciando correre le wl (nei range stabiliti) e poi lasciando correre la slit ma mediando le wl 


import numpy as np
import os, fnmatch
import matplotlib.pyplot as plt
import astropy
import idlsave
import scipy.io as spio
import matplotlib.ticker as ticker
from scipy import stats
import pylab


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




#all averaged slit
Msqy=np.mean(s.sq, axis=0)
Msuy=np.mean(s.su, axis=0)    #Msqy.shape (1240,) mediato su tutta la slit

#Blue wing I

polA0list= np.array([0,180,90])

if Msqy[a1] >0 and Msuy[a1] >=0:   
	polA0=polA0list[0]
elif Msqy[a1] >0 and Msuy[a1] <0:
	polA0=polA0list[1]
elif Msqy[a1] <0:
	polA0=polA0list[2] 
polA=0.5*((np.arctan2(Msuy[a1],Msqy[a1])* 180 / np.pi)+polA0) #arctan2(y, x) return arctan(y / x), in radians, "* 180 / np.pi" in degree.


for x in range(a1+1,b1+1): 
	if Msqy[x] >0 and Msuy[x] >=0:
		polA0=polA0list[0]
	elif Msqy[x] >0 and Msuy[x] <0:
		polA0=polA0list[1]
	elif Msqy[x] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msuy[x],Msqy[x])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])

plt.subplot(3,1,1)
ax = plt.gca()
plt.subplots_adjust(left  = 0.11, bottom=None, right = 0.92, top = 0.94, wspace=None, hspace = 0.5)
ax.scatter(wlarray[a1:b1+1],polAok,label='Blue wing I')
parafit=np.polyfit(wlarray[a1:b1+1],polAok,1) 
#print(parafit)
linfit=np.poly1d(parafit) 
ax.plot(wlarray[a1:b1+1],linfit(wlarray[a1:b1+1]), color='black')
ax.legend(loc='best')
ax.set_title('4"') 
ax.set_ylabel(r'$\alpha \ [degree]$')
ax.set_xlabel(r'$Wavelength\ [\AA]$')
ax.get_xaxis().get_major_formatter().set_useOffset(False)
	

#Red wing 

polA0list= np.array([0,180,90])

if Msqy[a2] >0 and Msuy[a2] >=0:
	polA0=polA0list[0]
elif Msqy[a2] >0 and Msuy[a2] <0:
	polA0=polA0list[1]
elif Msqy[a2] <0:
	polA0=polA0list[2] 
polA=0.5*((np.arctan2(Msuy[a2],Msqy[a2])* 180 / np.pi)+polA0) 


for x in range(a2+1,b2+1):  
	if Msqy[x] >0 and Msuy[x] >=0:
		polA0=polA0list[0]
	elif Msqy[x] >0 and Msuy[x] <0:
		polA0=polA0list[1]
	elif Msqy[x] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msuy[x],Msqy[x])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])


plt.subplot(3,1,2)
ax = plt.gca()
parafit=np.polyfit(wlarray[a2:b2+1],polAok,1) 
#print(parafit)
linfit=np.poly1d(parafit) 
plt.plot(wlarray[a2:b2+1],linfit(wlarray[a2:b2+1]), color='black')
ax.scatter(wlarray[a2:b2+1],polAok,label='Red wing')
ax.set_ylabel(r'$\alpha \ [degree]$')
ax.set_xlabel(r'$Wavelength\ [\AA]$')
ax.legend(loc='best') 
ax.get_xaxis().get_major_formatter().set_useOffset(False)
	

#Core 

polA0list= np.array([0,180,90])

if Msqy[ac] >0 and Msuy[ac] >=0:
	polA0=polA0list[0]
elif Msqy[ac] >0 and Msuy[ac] <0:
	polA0=polA0list[1]
elif Msqy[ac] <0:
	polA0=polA0list[2] 
polA=0.5*((np.arctan2(Msuy[ac],Msqy[ac])* 180 / np.pi)+polA0) 


for x in range(ac+1,bc+1):  
	if Msqy[x] >0 and Msuy[x] >=0:
		polA0=polA0list[0]
	elif Msqy[x] >0 and Msuy[x] <0:
		polA0=polA0list[1]
	elif Msqy[x] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msuy[x],Msqy[x])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])


plt.subplot(3,1,3)
ax = plt.gca()
parafit=np.polyfit(wlarray[ac:bc+1],polAok,1) 
#print(parafit)
linfit=np.poly1d(parafit) 
plt.plot(wlarray[ac:bc+1],linfit(wlarray[ac:bc+1]), color='black')
ax.scatter(wlarray[ac:bc+1],polAok,label='Core')
ax.set_ylabel(r'$\alpha \ [degree]$')
ax.set_xlabel(r'$Wavelength\ [\AA]$')
ax.legend(loc='best') 
ax.get_xaxis().get_major_formatter().set_useOffset(False)
	
plt.savefig('4227_m1_polA_all_av_slit.png')
plt.show() 

							

##################################
#polarization angle calcolato mediando lungo le wl e lasciando correre lungo la slit


#Blue wing I
sqx1=s.sq[0:140, a1:b1+1]   
Msqx1=np.mean(sqx1, axis=1) 

sux1=s.su[0:140, a1:b1+1] 
Msux1=np.mean(sux1, axis=1)    

#Red wing 
sqx2=s.sq[0:140, a2:b2+1]   
Msqx2=np.mean(sqx2, axis=1) 

sux2=s.su[0:140, a2:b2+1] 
Msux2=np.mean(sux2, axis=1)    

#Core 
sqxc=s.sq[0:140, ac:bc+1]   
Msqxc=np.mean(sqxc, axis=1) 

suxc=s.su[0:140, ac:bc+1] 
Msuxc=np.mean(suxc, axis=1)    

#Blue wing II :
sqx11=s.sq[0:140, a11:b11+1]   
Msqx11=np.mean(sqx11, axis=1) 

sux11=s.su[0:140, a11:b11+1] 
Msux11=np.mean(sux11, axis=1)    

#Blue wing III : 
sqx111=s.sq[0:140, a111:b111+1]   
Msqx111=np.mean(sqx111, axis=1) 

sux111=s.su[0:140, a111:b111+1] 
Msux111=np.mean(sux111, axis=1)    
  


#Blue wing I ....................................

polA0list= np.array([0,180,90])

if Msqx1[0] >0 and Msux1[0] >=0:
	polA0=polA0list[0]
elif Msqx1[0] >0 and Msux1[0] <0:
	polA0=polA0list[1]
elif Msqx1[0] <0:
	polA0=polA0list[2] 
polA=0.5*(np.arctan2(Msux1[0],Msqx1[0])* 180 / np.pi)+polA0


for y in range(1,140): 
	if Msqx1[y] >0 and Msux1[y] >=0:
		polA0=polA0list[0]
	elif Msqx1[y] >0 and Msux1[y] <0:
		polA0=polA0list[1]
	elif Msqx1[y] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msux1[y],Msqx1[y])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])

ax = plt.gca()
yrange=np.array(range(0,140))	
ax.plot(yrange,polAok,label='Blue wing I', color='blue')
ax.set_ylabel(r'$\alpha\ [degree]$')
ax.set_xlabel(r'$Slit\ position \ [pixel]$')
ax.set_ylim(-50,50)
ax.legend(loc='best')
ax.set_title('4"') 


#Red wing .....................................................

polA0list= np.array([0,180,90])

if Msqx2[0] >0 and Msux2[0] >=0:
	polA0=polA0list[0]
elif Msqx2[0] >0 and Msux2[0] <0:
	polA0=polA0list[1]
elif Msqx2[0] <0:
	polA0=polA0list[2] 
polA=0.5*(np.arctan2(Msux2[0],Msqx2[0])* 180 / np.pi)+polA0


for y in range(1,140): 
	if Msqx2[y] >0 and Msux2[y] >=0:
		polA0=polA0list[0]
	elif Msqx2[y] >0 and Msux2[y] <0:
		polA0=polA0list[1]
	elif Msqx2[y] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msux2[y],Msqx2[y])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])

ax.plot(yrange,polAok,label='Red wing', color='red')
ax.legend(loc='best') 


#Core.........................................................

polA0list= np.array([0,180,90])

if Msqxc[0] >0 and Msuxc[0] >=0:
	polA0=polA0list[0]
elif Msqxc[0] >0 and Msuxc[0] <0:
	polA0=polA0list[1]
elif Msqxc[0] <0:
	polA0=polA0list[2] 
polA=0.5*(np.arctan2(Msuxc[0],Msqxc[0])* 180 / np.pi)+polA0


for y in range(1,140): 
	if Msqxc[y] >0 and Msuxc[y] >=0:
		polA0=polA0list[0]
	elif Msqxc[y] >0 and Msuxc[y] <0:
		polA0=polA0list[1]
	elif Msqxc[y] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msuxc[y],Msqxc[y])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])

ax.plot(yrange,polAok,label='Core', color='black')
ax.legend(loc='best') 
	
plt.savefig('4227_m1_polA_wl_averaged.png')
plt.show()								

####################################################
#confronto solo tra blu wing I e red wing


#Blue wing I ....................................

polA0list= np.array([0,180,90])

if Msqx1[0] >0 and Msux1[0] >=0:
	polA0=polA0list[0]
elif Msqx1[0] >0 and Msux1[0] <0:
	polA0=polA0list[1]
elif Msqx1[0] <0:
	polA0=polA0list[2] 
polA=0.5*(np.arctan2(Msux1[0],Msqx1[0])* 180 / np.pi)+polA0


for y in range(1,140): 
	if Msqx1[y] >0 and Msux1[y] >=0:
		polA0=polA0list[0]
	elif Msqx1[y] >0 and Msux1[y] <0:
		polA0=polA0list[1]
	elif Msqx1[y] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msux1[y],Msqx1[y])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])

ax = plt.gca()
yrange=np.array(range(0,140))	
ax.plot(yrange,polAok,label='Blue wing I', color='blue')
ax.set_ylabel(r'$\alpha\ [degree]$')
ax.set_xlabel(r'$Slit\ position \ [pixel]$')
ax.set_ylim(-50,50)
ax.legend(loc='best') 
ax.set_title('4"')


#Red wing .....................................................

polA0list= np.array([0,180,90])

if Msqx2[0] >0 and Msux2[0] >=0:
	polA0=polA0list[0]
elif Msqx2[0] >0 and Msux2[0] <0:
	polA0=polA0list[1]
elif Msqx2[0] <0:
	polA0=polA0list[2] 
polA=0.5*(np.arctan2(Msux2[0],Msqx2[0])* 180 / np.pi)+polA0


for y in range(1,140): 
	if Msqx2[y] >0 and Msux2[y] >=0:
		polA0=polA0list[0]
	elif Msqx2[y] >0 and Msux2[y] <0:
		polA0=polA0list[1]
	elif Msqx2[y] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msux2[y],Msqx2[y])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])

ax.plot(yrange,polAok,label='Red wing', color='red')
ax.legend(loc='best') 



plt.savefig('4227_m1_polA_wl_averaged2.png')
plt.show()


######################################################################################################3
#--------------------confronto tra blue , blueI , blueII ----------------------------------

#Blue wing I ....................................

polA0list= np.array([0,180,90])

if Msqx1[0] >0 and Msux1[0] >=0:
	polA0=polA0list[0]
elif Msqx1[0] >0 and Msux1[0] <0:
	polA0=polA0list[1]
elif Msqx1[0] <0:
	polA0=polA0list[2] 
polA=0.5*(np.arctan2(Msux1[0],Msqx1[0])* 180 / np.pi)+polA0


for y in range(1,140): 
	if Msqx1[y] >0 and Msux1[y] >=0:
		polA0=polA0list[0]
	elif Msqx1[y] >0 and Msux1[y] <0:
		polA0=polA0list[1]
	elif Msqx1[y] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msux1[y],Msqx1[y])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])

ax = plt.gca()
yrange=np.array(range(0,140))	
ax.plot(yrange,polAok,label='Blue wing I', color='blue')
ax.set_ylabel(r'$\alpha [degree]$')
ax.set_xlabel(r'$Slit\ position\ [pixel]$')
ax.set_ylim(-50,50)
ax.legend(loc='best') 
ax.set_title('4"')

#Blue wing II ....................................

polA0list= np.array([0,180,90])

if Msqx11[0] >0 and Msux11[0] >=0:
	polA0=polA0list[0]
elif Msqx11[0] >0 and Msux11[0] <0:
	polA0=polA0list[1]
elif Msqx11[0] <0:
	polA0=polA0list[2] 
polA=0.5*(np.arctan2(Msux11[0],Msqx11[0])* 180 / np.pi)+polA0


for y in range(1,140): 
	if Msqx11[y] >0 and Msux11[y] >=0:
		polA0=polA0list[0]
	elif Msqx11[y] >0 and Msux11[y] <0:
		polA0=polA0list[1]
	elif Msqx11[y] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msux11[y],Msqx11[y])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])

ax.plot(yrange,polAok,label='Blue wing II', color='black')
#ax.set_ylim(-50,50)
ax.legend(loc='best') 


#Blue wing III ....................................

polA0list= np.array([0,180,90])

if Msqx111[0] >0 and Msux111[0] >=0:
	polA0=polA0list[0]
elif Msqx111[0] >0 and Msux111[0] <0:
	polA0=polA0list[1]
elif Msqx111[0] <0:
	polA0=polA0list[2] 
polA=0.5*(np.arctan2(Msux111[0],Msqx111[0])* 180 / np.pi)+polA0


for y in range(1,140): 
	if Msqx111[y] >0 and Msux111[y] >=0:
		polA0=polA0list[0]
	elif Msqx111[y] >0 and Msux111[y] <0:
		polA0=polA0list[1]
	elif Msqx111[y] <0:
		polA0=polA0list[2] 
	polA=np.append(polA,0.5*(np.arctan2(Msux111[y],Msqx111[y])* 180 / np.pi)+polA0)

#reimposto l'angolo di pol tra -90 e 90
polAok=polA[0]  
for j in range(1,len(polA)):
	if polA[j]>90:
		polAok=np.append(polAok, polA[j]-180)
	else:
		polAok=np.append(polAok, polA[j])

ax.plot(yrange,polAok,label='Blue wing III', color='green')
#ax.set_ylim(-50,50)
ax.legend(loc='best') 
plt.savefig('4227_m1_polA_wltris_averaged.png')
plt.show()								




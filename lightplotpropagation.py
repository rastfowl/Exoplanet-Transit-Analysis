# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:34:55 2018

@author: rastf
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Loading in Data - delta, epsilon, zeta and eta are comparison stars. Q is main star
# =============================================================================
data = np.loadtxt('rawdatasedgewick.csv', delimiter=',') 
number      =data[:,0] 
JD	        =data[:,3]
AIRMASS	    =data[:,10]
fluxQ       =data[:,16]
fluxdelta   =data[:,17]
fluxepsil   =data[:,18]
fluxzeta    =data[:,19]
fluxeta     =data[:,20]
fluxerrQ    =data[:,21]	
fluxerrdelta=data[:,22]	
fluxerrepsil=data[:,23]	
fluxerrzeta =data[:,24]
fluxerreta  =data[:,25]	
SNRQ        =data[:,26]	
SNRdelta    =data[:,27]	
SNRepsil    =data[:,28]	
SNRzeta     =data[:,29]	
SNReta      =data[:,30]	
# =============================================================================
# Creating dot-plot of Star Q, as well as the same with flux offset of comparisons
# =============================================================================
plt.plot(JD, fluxQ-fluxdelta, '.', label='Star Q - Star $\delta$')
plt.plot(JD, fluxQ-fluxepsil-0.17, '.', label='Star Q - Star $\epsilon$')
plt.plot(JD, fluxQ-fluxzeta-0.2, '.', label='Star Q - Star $\zeta$')
plt.plot(JD, fluxQ-fluxeta-0.27, '.', label='Star Q - Star $\eta$')
plt.legend(loc='best')
plt.xlabel('Time (Julian Date J.D.)')
plt.ylabel('Flux')
plt.title('Flux of Qatar-1, with comparison stars as offsets')
# =============================================================================
# Creating line plot of fluxes of Star Q and comparisons
# =============================================================================
plt.figure()
plt.plot(JD,fluxdelta-0.006,label='Star $\delta$')
plt.plot(JD,fluxepsil+0.165,label='Star $\epsilon$')
plt.plot(JD,fluxeta+0.278,label='Star $\eta$')
plt.plot(JD,fluxzeta+0.2,label='Star $\zeta$')
plt.plot(JD,fluxQ-3.5,label='Qatar-1')
plt.xlabel('Julian Date')
plt.ylabel('Flux (arbitrary units)')
plt.legend(loc='best')
# =============================================================================
# Creating an averaged-out light curve
# =============================================================================
testy=np.zeros(len(JD))
testy2=np.zeros(len(JD))
testy3=np.zeros(len(JD))
testy4=np.zeros(len(JD))
testy5=np.zeros(len(JD))
testy6=np.zeros(len(JD))
for i in range(len(JD)):
    testy[i]=(fluxQ[i-1]+fluxQ[i])*0.5
for i in range(len(JD)):
    testy2[i]=(testy[i-1]+testy[i])*0.5
for i in range(len(JD)):
    testy3[i]=(testy2[i-1]+testy2[i])*0.5
for i in range(len(JD)):
    testy4[i]=(testy3[i-1]+testy3[i])*0.5
for i in range(len(JD)):
    testy5[i]=(testy4[i-1]+testy4[i])*0.5
for i in range(len(JD)):
    testy6[i]=(testy5[i-1]+testy5[i])*0.5
plt.figure()
plt.plot(JD,fluxQ,label='orig')
plt.plot(JD,testy,label='1st')
plt.plot(JD,testy2,label='2nd')
plt.plot(JD,testy3,label='3rd')
plt.plot(JD,testy4,label='4th')
plt.plot(JD,testy5,label='5th')
plt.plot(JD,testy6,label='6th')
plt.legend(loc='best')
# =============================================================================
# FLUX ---> MAGNITUDE
# =============================================================================
magQ=-2.5*(np.log10(fluxQ))
magdelta=-2.5*(np.log10(fluxdelta))
magepsil=-2.5*(np.log10(fluxepsil))
magzeta=-2.5*(np.log10(fluxzeta))
mageta=-2.5*(np.log10(fluxeta))
magQerr=2.5*(fluxerrQ)/fluxQ
# =============================================================================
# CREATING ARRAY TO BE PUT INTO ETD
# =============================================================================
t=np.zeros((181,3))
t[:,0]=JD
t[:,1]=magQ
t[:,2]=magQerr
np.savetxt('putmeinetd.txt', t)
# =============================================================================
# LOADING IN ETD FIT
# =============================================================================
data2=np.loadtxt('fitdatasedgewick.txt')
JD2=data2[:,0] # JD == JD2 + const (55796)
magQ2=data2[:,1] # mag2 == magQ
magfit=data2[:,2] 
#plotting the transit fit against the data
def flux(mag):
    return 10**(-mag/2.5)
fluxfit=flux(magfit)
# =============================================================================
# Plotting the ETD fit
# =============================================================================
plt.figure()
plt.errorbar(JD2, fluxQ, fmt='g.', yerr=fluxerrQ, ecolor='green', capsize=3, errorevery=1, elinewidth=1, label='Flux of Qatar-1')
plt.plot(JD2,fluxfit, label='Trend fit', linewidth=3)
plt.ylabel('Flux')
plt.xlabel('Julian Date')
plt.title('Light curve for transit of Qatar-1b')
plt.legend(loc='best')
plt.grid()
plt.figure()
# =============================================================================
# Finding first and second derivatives
# =============================================================================
grads=np.zeros(len(JD))
for i in range(len(JD)):
    grads[i]=(magfit[i]-magfit[i-1])/(JD2[i]-JD2[i-1])
#calculating the second derivative
secderiv=np.zeros(len(JD))
for i in range(len(JD)):
    secderiv[i]=(grads[i]-grads[i-1])/(JD2[i]-JD2[i-1])
#calculating scaling factors so the first and second derivatives of the fluxfit can be plotted on top of the fluxfit
avgflux=sum(fluxfit)/len(fluxfit)
heightgrads=(max(fluxfit)-min(fluxfit))/(max(grads)-min(grads))
heightsecderiv=(max(fluxfit)-min(fluxfit))/(max(secderiv)-min(secderiv))
#plotting the fluxfit against its derivatives
plt.plot(JD2,flux(magfit), label='ETD Light-curve fit', linewidth=3)
plt.plot(JD2,heightgrads*grads+avgflux, label='First derivative of fitting')
plt.plot(JD2,heightsecderiv*secderiv+avgflux, label='Second derivative of fitting')
plt.ylabel('Flux')
plt.xlabel('Julian Date')
plt.legend(loc='best')
# =============================================================================
# Finding times of beginning and ending of ingress and egress
# =============================================================================
def timewhen(func,thing):
    for i in range(len(func)):
        if func[i] == thing:
            return JD[i]
def indexwhen(func,thing):
    for i in range(len(func)):
        if func[i] == thing:
            return i
# mid == middle of transit
# time1 == beginning of ingress
# time2 == end of ingress
# time3 == beginning of egress
# time4 == end of egress
# depth flux == depth of curve in flux
# depth mag == depth of curve in magnitudes
# FINDING MID
mid = timewhen(fluxfit,min(fluxfit))
midindex = indexwhen(fluxfit,min(fluxfit))
#Splitting 1stder and 2ndder into two parts
gradsfirsthalf=grads[:midindex]
gradssecndhalf=grads[midindex:]
secderivfirsthalf=secderiv[:midindex]
secderivsecndhalf=secderiv[midindex:]

# FINDING TIME1
time1=timewhen(secderiv,max(secderivfirsthalf))
time1index=indexwhen(secderiv,max(secderivfirsthalf))
time1=time1*86400
# FINDING TIME2
time2=timewhen(secderiv,min(secderivfirsthalf))
time2index=indexwhen(secderiv,min(secderivfirsthalf))
time2=time2*86400
# FINDING TIME3
time3=timewhen(secderiv,min(secderivsecndhalf))
time3index=indexwhen(secderiv,min(secderivsecndhalf))
time3=time3*86400
# FINDING TIME4
time4=timewhen(secderiv,max(secderivsecndhalf))
time4index=indexwhen(secderiv,max(secderivsecndhalf))
time4=time4*86400

# =============================================================================
# Radius of planet
# =============================================================================
# FINDING DEPTHFLUX
fluxpolyfit=np.polyfit(JD2[time4index:],fluxfit[time4index:],deg=1)
bestfitflux=fluxpolyfit[0]*JD2+fluxpolyfit[1]
depthflux=1-(fluxfit[midindex]/bestfitflux[midindex])
fluxerr=fluxerrQ[midindex]
fluxmin=fluxfit[midindex]
fluxmax=bestfitflux[midindex]
solarradius=695508000
qatar1radius=(0.803)*solarradius
qatar1radiuserr=0.01*solarradius
radjup=69911000
radjuperr=1

fminterm=((qatar1radius)*((1-fluxmin/fluxmax)**(-0.5))*(1/(2*fluxmax)))**2*(fluxerr**2)
fmaxterm=((qatar1radius)*((1-fluxmin/fluxmax)**(-0.5))*(1/(2*(fluxmax**2))))**2*(fluxerr**2)
Rqterm=(1)*((1-fluxmin/fluxmax)**0.5)*(qatar1radiuserr**2)

radiuserr=(fminterm+fmaxterm+Rqterm)**0.5
radius=((1-(fluxmin/fluxmax))**0.5)*qatar1radius

radiuserrjup=radiuserr/radjup
radiusjup=radius/radjup


# =============================================================================
# Orbital speed of planet
# =============================================================================
timeerr=-(sum(JD2[:-1]-JD2[1:])/len(JD2))*86400
radterm=((2/(time2-time1))**2)*((radiuserr)**2)
timeterm=2*((2*radius/((time2-time1)**2))**2)*(timeerr**2)
velocityerr=(radterm+timeterm)**0.5
velocity=2*radius/((time2-time1))
velocitykm=velocity/1000
velocitykmerr=velocityerr/1000

stellarspeedkm=0.218
stellarspeedkmerr=0.016


# =============================================================================
# Mass of planet
# =============================================================================
solarmass=2*(10**30)
qatar1mass=0.838*solarmass
qatar1masserr=0.016*solarmass

massterm=(qatar1masserr**2)*(stellarspeedkm/velocitykm)**2
stelvelterm=(stellarspeedkmerr**2)*(qatar1mass/velocitykm)**2
planvelterm=(velocitykmerr**2)*(qatar1mass*stellarspeedkm/(velocitykm**2))**2
masserr=(massterm+stelvelterm+planvelterm)**0.5


mass=qatar1mass*stellarspeedkm/velocitykm

jup=1.898*(10**27)
massjup=mass/jup
massjuperr=masserr/jup

# =============================================================================
# Density of planet
# =============================================================================
density=0.001*0.75*(np.pi**-1)*mass*(radius**-3)

densmassterm=((0.75*(np.pi**-1)*(radius**-3))**2)*(masserr**2)
densradterm=((2.25*(np.pi**-1)*(radius**-4)*mass)**2)*(radiuserr**2)
densityerr=(densmassterm+densradterm)**0.5

planetmasses=np.array([0.330,4.87,5.97,0.642,1898,568,86.8,102,(mass/(10**24))])*(10**24)
planetdiamscaled=np.array([4879,12104,12756,6792,142984,120536,51118,49528,((2*radius)/1000)])/10000
planetdensities=np.array([5.427,5.243,5.514,3.933,1.326,0.687,1.271,1.638,density])

plt.figure()
#plt.plot(planetmasses[0],planetdensities[0], '.', label='Mercury', markersize=int(20*planetdiamscaled[0]))
#plt.plot(planetmasses[1],planetdensities[1], '.', label='Venus', markersize=int(20*planetdiamscaled[1]))
#plt.plot(planetmasses[2],planetdensities[2], '.', label='Earth', markersize=int(20*planetdiamscaled[2]))
#plt.plot(planetmasses[3],planetdensities[3], '.', label='Mars', markersize=int(20*planetdiamscaled[3]))
plt.plot(planetmasses[4],planetdensities[4], 'r.', label='Jupiter', markersize=int(7*planetdiamscaled[4]))
plt.plot(planetmasses[5],planetdensities[5], 'm.', label='Saturn', markersize=int(7*planetdiamscaled[5]))
plt.plot(planetmasses[6],planetdensities[6], 'c.', label='Uranus', markersize=int(7*planetdiamscaled[6]))
plt.plot(planetmasses[7],planetdensities[7], 'b.', label='Neptune', markersize=int(7*planetdiamscaled[7]))
plt.plot(planetmasses[8],planetdensities[8], 'k.', label='Qatar-1b', markersize=int(7*planetdiamscaled[8]))
plt.xlabel('Mass ($10^{28} kg$)')
plt.ylabel('Density ($g/cm^3$)')
plt.xlim(right=2.1*10**27)
plt.ylim(bottom=0.5)
plt.legend(loc='best', markerscale=0.1)


print('Radius is %.3e +/- %.3e m.\n    This is equivalent to %.2f +/- %.2f R(Jup).' % (radius,radiuserr,radiusjup,radiuserrjup))
print('Mass is %.3e +/- %.3e kg.\n    This is equivalent to %.2f +/- %.2f M(Jup).' % (mass,masserr,massjup,massjuperr))





















#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 16:41:49 2017

@author: xuduo
"""





import numpy as np
import matplotlib.pyplot as plt
import astropy
import astropy.constants









def density_profile_single_pl(M_disk,p,distance):
    R_in=6*AU
    R_out=310*AU
    
    surf_den_100=M_disk*(2-p)/(2*np.pi*(AU*100)**p*(R_out**(2-p)-R_in**(2-p)))
#    surf_den_100=M_disk*(2-p)/(2*np.pi(100)**p*(R_out**(2-p)-R_in**(2-p)))
    
    surf_den=surf_den_100*(distance/(100*AU))**(-p)
#    surf_den=surf_den_100*(distance/(100))**(-p)

    return surf_den


def density_profile_double_pl(M_disk,distance):
    R_in=6.7*AU
    R_out=310*AU
    R_T=88*AU
    p1=-2.7
    p2=1.5
    m1=R_T**(2-p1)-R_in**(2-p1)
    m2=R_out**(2-p2)-R_T**(2-p2)
    
    surf_den_100=M_disk*(2-p1)*(2-p2)/(2*np.pi*(m1*(2-p2)*R_T**p1+m2*(2-p1)*R_T**p2))
#    surf_den_100=M_disk*(2-p)/(2*np.pi(100)**p*(R_out**(2-p)-R_in**(2-p)))
    
    surf_den=surf_den_100*(distance/(100*AU))**(-p)
#    surf_den=surf_den_100*(distance/(100))**(-p)

    return surf_den


def surface_den_to_number(surf_den,R_dust):
    surf_num_den=surf_den/(4./3.*np.pi*rho_dust*R_dust**3)
    return surf_num_den



R_sun=astropy.constants.R_sun.cgs.to_value()
L_sun=astropy.constants.L_sun.cgs.to_value()
M_sun=astropy.constants.M_sun.cgs.to_value()
sigma_sb=astropy.constants.sigma_sb.cgs.to_value()
AU=astropy.constants.au.cgs.to_value()
pc=astropy.constants.pc.cgs.to_value()
rho_dust=2    

R_dust_01=0.1*1e-4
R_dust_1=1*1e-4
R_dust_10=10*1e-4
R_dust_1m=0.1
    
star_radius=R_sun*1.88
star_temp=8790#8590
M_disk=10**(-1.16) * M_sun
p=1.29
#p=1.5
distance=np.linspace(6,310,200)  * AU

###constant
c=2.99792458*10**10
k=1.380650*10**(-16)
h=6.626069*10**(-27)
e=np.exp(1)

###variable 
v=np.logspace(11.4,15.8,num=10000)
#v=np.logspace(10,16,num=10000)
#v_ghz=np.logspace(9,13,num=1000)/1e9
lambda_1=c/v*1e4


surf_density_0=density_profile_single_pl(M_disk,1,distance)
surf_density_1=density_profile_single_pl(M_disk,p,distance)
surf_density_2=density_profile_single_pl(M_disk,1.5,distance)
surf_density_3=density_profile_single_pl(M_disk,1.9,distance)
surf_density_4=density_profile_double_pl(M_disk,distance)


plt.figure(1)
plt.clf()
plt.plot(distance/AU,surf_density_0,label='1')
plt.plot(distance/AU,surf_density_1,label='1.29')
plt.plot(distance/AU,surf_density_2,label='1.5')
plt.plot(distance/AU,surf_density_3,label='1.9')
plt.plot(distance/AU,surf_density_4,label='dpl')
plt.legend()
plt.xlabel(r'$\rm distance\, (AU)$')
plt.ylabel(r'$\rm surface\, density\, (g/cm^{-2})$')
plt.yscale('log')

plt.savefig('surface_density.pdf',bbox_inches='tight')




plt.figure(2)
plt.clf()
plt.plot(distance/AU,surface_den_to_number(surf_density_0,R_dust_01),label='1')
plt.plot(distance/AU,surf_density_1,label='1.29')
plt.plot(distance/AU,surf_density_2,label='1.5')
plt.plot(distance/AU,surf_density_3,label='1.9')
plt.plot(distance/AU,surf_density_4,label='dpl')
plt.legend()
plt.xlabel(r'$\rm distance\, (AU)$')
plt.ylabel(r'$\rm surface\, number \, density\, (cm)$')
plt.yscale('log')

#plt.savefig('surface_number_density.pdf',bbox_inches='tight')





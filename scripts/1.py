#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 02:06:54 2017

@author: xuduo
"""


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy
import astropy.constants
from scipy.integrate import simps
from scipy.interpolate import interp1d
from lmfit import minimize, Parameters
import lmfit


###plank function
def brightness_temp(temp,v):
    b_t=2*h*v**3/c**2/(e**(h*v/k/temp)-1)#*1e23  ## in Jy
#    b_t=2*h*c**2/lambda_1**5/(e**(h*c/k/temp/lambda_1)-1)*1e23
    return b_t

### dust model
def dust_model_Qabs(dust_model_01,lambda_1):
    lambda_dust_01=dust_model_01[:,0]
    Qabs_01=dust_model_01[:,1]
    f2 = interp1d(lambda_dust_01,Qabs_01, kind='slinear',fill_value='extrapolate')   ### 'cubic'
    Qabs_01_inter=f2(lambda_1) 
    return Qabs_01_inter

def dust_power_out(Qabs_01_inter,temp_low,temp_high,v,R_dust_01,P_in):
    
#    P_out_1=np.pi*simps(Qabs_01_inter*brightness_temp(temp_low,v)*np.pi*R_dust_01**2,v,dx=1)
#    P_out_2=np.pi*simps(Qabs_01_inter*brightness_temp(temp_high,v)*np.pi*R_dust_01**2,v,dx=1)
    temp_middle=(temp_low+temp_high)/2.0
    P_out_middle=np.pi*simps(Qabs_01_inter*brightness_temp(temp_middle,v)*1e23*np.pi*4*R_dust_01**2,v,dx=1)
    if P_out_middle<P_in:
        return temp_middle,temp_high
    if P_out_middle>P_in:
        return temp_low,temp_middle


def dust_temp_fit(Qabs_01_inter,temp_low1,temp_high1,v,R_dust_01,P_in_01,delta_num,ctt_num):
    temp_low=temp_low1
    temp_high=temp_high1
    delta=temp_high-temp_low
    ctt=0
    while (delta>delta_num) and (ctt<ctt_num):
        temp_low_new,temp_high_new=dust_power_out(Qabs_01_inter,temp_low,temp_high,v,R_dust_01,P_in_01)
        temp_low=temp_low_new
        temp_high=temp_high_new
        delta=temp_high-temp_low
        ctt=ctt+1
#    return temp_low_new,temp_high_new
    return (temp_low+temp_high)/2.0

    
def dust_spec(dust_temp,Qabs_01_inter,v,R_dust_01):
    b_t_dust=brightness_temp(dust_temp,v)*1e23
    L_nu_dust=np.pi*b_t_dust*4*np.pi*R_dust_01**2
    f_nu_dust_1=L_nu_dust*Qabs_01_inter
#    f_nu_dust_1=f_nu_dust_1/(4./3*np.pi*R_dust_01**3*rho_dust)
    return f_nu_dust_1
    
def convole_band_model(dust_spec_1,filter_data):
    f2 = interp1d(lambda_1,dust_spec_1, kind='slinear',fill_value='extrapolate')   ### 'cubic'
    lambda_obs=filter_data[:,0]/1e4
    spec_01_1_inter=f2(lambda_obs) 
    intensity_band=simps(spec_01_1_inter*filter_data[:,1],lambda_obs,dx=0.1)
    return intensity_band

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


def append_data(flux_4,flux_5,flux_8,flux_24,flux_70,flux_160,flux_870,dust_spec_01_1,R_dust):
    flux_4.append(convole_band_model(dust_spec_01_1,fliter_4)*surface_den_to_number(surf_density_1,R_dust))
    flux_5.append(convole_band_model(dust_spec_01_1,fliter_5)*surface_den_to_number(surf_density_1,R_dust))
    flux_8.append(convole_band_model(dust_spec_01_1,fliter_8)*surface_den_to_number(surf_density_1,R_dust))
    flux_24.append(convole_band_model(dust_spec_01_1,fliter_24)*surface_den_to_number(surf_density_1,R_dust))
    flux_70.append(convole_band_model(dust_spec_01_1,fliter_70)*surface_den_to_number(surf_density_1,R_dust))
    flux_160.append(convole_band_model(dust_spec_01_1,fliter_160)*surface_den_to_number(surf_density_1,R_dust))
    flux_870.append(convole_band_model(dust_spec_01_1,fliter_870)*surface_den_to_number(surf_density_1,R_dust))
    return flux_4,flux_5,flux_8,flux_24,flux_70,flux_160,flux_870
    
def append_data_2(flux_4_01,flux_5_01,flux_8_01,flux_24_01,flux_70_01,flux_160_01,flux_870_01):
    new_array=np.zeros([200,7])
    flux_4_01=np.asarray(flux_4_01)
    flux_5_01=np.asarray(flux_5_01)
    flux_8_01=np.asarray(flux_8_01)
    flux_24_01=np.asarray(flux_24_01)
    flux_70_01=np.asarray(flux_70_01)
    flux_160_01=np.asarray(flux_160_01)
    flux_870_01=np.asarray(flux_870_01)
    new_array[:,0]=flux_4_01
    new_array[:,1]=flux_5_01
    new_array[:,2]=flux_8_01
    new_array[:,3]=flux_24_01
    new_array[:,4]=flux_70_01
    new_array[:,5]=flux_160_01
    new_array[:,6]=flux_870_01
    
    return    new_array 
    
    
R_sun=astropy.constants.R_sun.cgs.to_value()
L_sun=astropy.constants.L_sun.cgs.to_value()
M_sun=astropy.constants.M_sun.cgs.to_value()
sigma_sb=astropy.constants.sigma_sb.cgs.to_value()
AU=astropy.constants.au.cgs.to_value()
pc=astropy.constants.pc.cgs.to_value()
rho_dust=2    


star_radius=R_sun*1.88
star_temp=8790#8590
M_disk=10**(-1.16) * M_sun
#p=1.29
p=1.9
distance=np.linspace(6,310,200)  * AU


fliter_4=np.loadtxt('../filter/Spitzer-IRAC.I2.dat')
fliter_5=np.loadtxt('../filter/Spitzer-IRAC.I3.dat')
fliter_8=np.loadtxt('../filter/Spitzer-IRAC.I4.dat')
fliter_24=np.loadtxt('../filter/Spitzer-MIPS.24mu.dat')
fliter_70=np.loadtxt('../filter/Spitzer-MIPS.70mu.dat')
fliter_160=np.loadtxt('../filter/Spitzer-MIPS.160mu.dat')
fliter_870=np.zeros((10,2))
fliter_870[:,0]=np.linspace(862.8,883,10)*1e4
fliter_870[:,1]=1


L_star=4*np.pi*star_radius**2*sigma_sb*star_temp**4


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



#"""
dust_model_01=np.loadtxt('../dust_model/0.1um.txt')
dust_model_1=np.loadtxt('../dust_model/1um.txt')
dust_model_10=np.loadtxt('../dust_model/10um.txt')


#surf_density=density_profile_single_pl(M_disk,p,distance)
surf_density=density_profile_double_pl(M_disk,distance)


#plt.figure(1)
#plt.clf()
#plt.plot(distance/AU,surf_density)
#plt.xlabel(r'$\rm distance\, (AU)$')
#plt.ylabel(r'$\rm surface density\, (g\, cm^{-2})$')

flux_4_01=[]
flux_5_01=[]
flux_8_01=[]
flux_24_01=[]
flux_70_01=[]
flux_160_01=[]
flux_870_01=[]

flux_4_1=[]
flux_5_1=[]
flux_8_1=[]
flux_24_1=[]
flux_70_1=[]
flux_160_1=[]
flux_870_1=[]

flux_4_10=[]
flux_5_10=[]
flux_8_10=[]
flux_24_10=[]
flux_70_10=[]
flux_160_10=[]
flux_870_10=[]

flux_4_1m=[]
flux_5_1m=[]
flux_8_1m=[]
flux_24_1m=[]
flux_70_1m=[]
flux_160_1m=[]
flux_870_1m=[]

temp_all=np.zeros((200,4))

for ctt_1,Distance_1,surf_density_1 in zip(range(200),distance,surf_density):
    
    Qabs_01_inter=dust_model_Qabs(dust_model_01,lambda_1)
    Qabs_1_inter=dust_model_Qabs(dust_model_1,lambda_1)
    Qabs_10_inter=dust_model_Qabs(dust_model_10,lambda_1)
    Qabs_1m=np.ones(len(v))
    R_dust_01=0.1*1e-4
    R_dust_1=1*1e-4
    R_dust_10=10*1e-4
    R_dust_1m=0.1
    
    
    b_t_star=brightness_temp(star_temp,v)
    L_nu_star=np.pi*b_t_star*4*np.pi*star_radius**2
    f_nu_star_1=L_nu_star/(4*np.pi*(Distance_1)**2)*1e23
    
    P_in_01=simps(Qabs_01_inter*f_nu_star_1*np.pi*R_dust_01**2,v,dx=1)
    P_in_1=simps(Qabs_1_inter*f_nu_star_1*np.pi*R_dust_1**2,v,dx=1)
    P_in_10=simps(Qabs_10_inter*f_nu_star_1*np.pi*R_dust_10**2,v,dx=1)
    P_in_1m=simps(Qabs_1m*f_nu_star_1*np.pi*R_dust_1m**2,v,dx=1)
    
    
    T_01_1= dust_temp_fit(Qabs_01_inter,1,1e4,v,R_dust_01,P_in_01,0.05,20)
    T_1_1= dust_temp_fit(Qabs_1_inter,1,1e4,v,R_dust_1,P_in_1,0.05,20)
    T_10_1= dust_temp_fit(Qabs_10_inter,1,1e4,v,R_dust_10,P_in_10,0.05,20)
    T_1m_1= dust_temp_fit(Qabs_1m,1,1e4,v,R_dust_1m,P_in_1m,0.05,20)
    
    temp_all[ctt_1,0]=T_01_1
    temp_all[ctt_1,1]=T_1_1
    temp_all[ctt_1,2]=T_10_1
    temp_all[ctt_1,3]=T_1m_1

        
    dust_spec_01_1=dust_spec(T_01_1,Qabs_01_inter,v,R_dust_01)*1e23
    dust_spec_1_1=dust_spec(T_1_1,Qabs_1_inter,v,R_dust_1)*1e23
    dust_spec_10_1=dust_spec(T_10_1,Qabs_10_inter,v,R_dust_10)*1e23
    dust_spec_1m_1=dust_spec(T_1m_1,Qabs_1m,v,R_dust_1m)*1e23
    
    
    flux_4_01,flux_5_01,flux_8_01,flux_24_01,flux_70_01,flux_160_01,flux_870_01=append_data(flux_4_01,flux_5_01,flux_8_01,flux_24_01,flux_70_01,flux_160_01,flux_870_01,dust_spec_01_1,R_dust_01)    
    flux_4_1,flux_5_1,flux_8_1,flux_24_1,flux_70_1,flux_160_1,flux_870_1=append_data(flux_4_1,flux_5_1,flux_8_1,flux_24_1,flux_70_1,flux_160_1,flux_870_1,dust_spec_1_1,R_dust_1)    
    flux_4_10,flux_5_10,flux_8_10,flux_24_10,flux_70_10,flux_160_10,flux_870_10=append_data(flux_4_10,flux_5_10,flux_8_10,flux_24_10,flux_70_10,flux_160_10,flux_870_10,dust_spec_10_1,R_dust_10)    
    flux_4_1m,flux_5_1m,flux_8_1m,flux_24_1m,flux_70_1m,flux_160_1m,flux_870_1m=append_data(flux_4_1m,flux_5_1m,flux_8_1m,flux_24_1m,flux_70_1m,flux_160_1m,flux_870_1m,dust_spec_1m_1,R_dust_1m)    
    



flux_all_01=append_data_2(flux_4_01,flux_5_01,flux_8_01,flux_24_01,flux_70_01,flux_160_01,flux_870_01)
flux_all_1=append_data_2(flux_4_1,flux_5_1,flux_8_1,flux_24_1,flux_70_1,flux_160_1,flux_870_1)
flux_all_10=append_data_2(flux_4_10,flux_5_10,flux_8_10,flux_24_10,flux_70_10,flux_160_10,flux_870_10)
flux_all_1m=append_data_2(flux_4_1m,flux_5_1m,flux_8_1m,flux_24_1m,flux_70_1m,flux_160_1m,flux_870_1m)

#plt.plot(distance,(flux_4)/np.median(flux_4))
#plt.plot(distance,(flux_5)/np.median(flux_5))
#plt.plot(distance,(flux_8)/np.median(flux_8))
#plt.plot(distance,(flux_24)/np.median(flux_24))
#plt.plot(distance,(flux_70)/np.median(flux_70))
#plt.plot(distance,(flux_160)/np.median(flux_160))
#plt.plot(distance,(flux_870)/np.median(flux_870))
#plt.plot(distance,flux_160_01,label='01')
#plt.plot(distance,flux_160_10,label='1')
#plt.plot(distance,flux_160_1,label='10')
#plt.plot(distance,flux_160_1m,label='1m')
#
#
#plt.plot(distance,flux_870_01,label='01')
#plt.plot(distance,flux_870_10,label='1')
#plt.plot(distance,flux_870_1,label='10')
#plt.plot(distance,flux_870_1m,label='1m')


new_array=np.zeros((200,7,4))
new_array[:,:,0]=flux_all_01
new_array[:,:,1]=flux_all_1
new_array[:,:,2]=flux_all_10
new_array[:,:,3]=flux_all_1m

#np.save('temp_p_'+str(p)+'_single_pw.npy',temp_all)
#np.save('p_'+str(p)+'_single_pw.npy',new_array)
np.save('temp_double_pw.npy',temp_all)
np.save('double_pw.npy',new_array)





fig = plt.figure(2)
fig.clf()
ax1 = fig.add_subplot(111)




xlim=np.array([3,1e3])
ylim=np.array([1e-20,1e36])

ax1.plot(lambda_1,dust_spec_01_1,label="dust 0.1um",color="green",linewidth=2)
ax1.plot(lambda_1,dust_spec_1_1,label="dust 1um",color="blue",linewidth=2)
ax1.plot(lambda_1,dust_spec_10_1,label="dust 10um",color="orange",linewidth=2)
ax1.plot(lambda_1,dust_spec_1m_1,label="dust 1mm",color="red",linewidth=2)



ax1.set_xlim(xlim)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.yaxis.get_data_interval()


ax2 = ax1.twiny()
ax2.set_xscale('log')
ax2.set_xlabel(r"$\nu\, [Hz]$",fontsize=15)
ax2.set_xlim(c/xlim*1e4)
ax2.set_ylim(ylim)

#ax2.set_xticks([1e11,1e10,1e9])
#ax2.set_xticklabels([r'$10^1$','8','99','1'])



#plt.xlim(9,4.5e3)
#plt.ylim(1e-19,1e-11)
#ax1.tick_params(labelsize=14)


ax1.set_xlabel(r"$\lambda\, [\mu m]$",fontsize=15)
ax1.set_ylabel(r'$\rm L_{\nu}\, [Jy\, cm^{2}]$',fontsize=15)
#ax1.set_title("Blackbody",fontsize=14)
ax1.legend(loc=0,fontsize=14)
#plt.savefig('part_3_dust_sed.pdf',bbox_inches='tight')
#fig.suptitle('310 AU')
plt.savefig('sed.pdf',bbox_inches='tight')
plt.show()
#"""



#
#for ctt_1 in range(4):
#    plt.plot(distance/AU,new_array[:,0,ctt_1])
#plt.title(str(p))
#plt.yscale('log')


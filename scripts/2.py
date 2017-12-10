#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 17:31:03 2017

@author: xuduo
"""


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import astropy
import astropy.constants
from scipy.integrate import simps
from scipy.interpolate import interp1d




def  generate_image(data_all):
    a, b = 310, 310
    n = 620
    new_array_all=np.zeros((620,620,7,4,5))    
    y,x = np.ogrid[-a:n-a, -b:n-b]
    
    for ctt in range(1,len(distance)):
        mask = (x*x + y*y >= distance[ctt-1]**2)&(x*x + y*y < distance[ctt]**2)
        for ctt1 in range(7):
            for ctt2 in range(4):
                for ctt3 in range(5):
                    new_array_all[mask,ctt1,ctt2,ctt3] = data_all[ctt,ctt1,ctt2,ctt3]
    return new_array_all



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



double_pw=np.load('double_pw.npy')
p_1_single_pw=np.load('p_1_single_pw.npy')
p_15_single_pw=np.load('p_1.5_single_pw.npy')
p_19_single_pw=np.load('p_1.9_single_pw.npy')
p_129_single_pw=np.load('p_1.29_single_pw.npy')


data_all=np.zeros((200,7,4,5))
data_all[:,:,:,0]=p_1_single_pw
data_all[:,:,:,1]=p_15_single_pw
data_all[:,:,:,2]=p_19_single_pw
data_all[:,:,:,3]=p_129_single_pw
data_all[:,:,:,4]=double_pw




distance=distance/AU
#a, b = 310, 310
#n = 620
#array_01 = np.zeros((n, n))
#array_1 = np.zeros((n, n))
#array_10 = np.zeros((n, n))
#array_1m = np.zeros((n, n))
#
#array_sp1 = np.zeros((n, n))
#array_sp2 = np.zeros((n, n))
#array_sp3 = np.zeros((n, n))
#array_sp4 = np.zeros((n, n))
#array_dp1 = np.zeros((n, n))
#
#y,x = np.ogrid[-a:n-a, -b:n-b]
#
#for ctt in range(1,len(distance)):
#    mask = (x*x + y*y >= distance[ctt-1]**2)&(x*x + y*y < distance[ctt]**2)
##    array_01[mask] = double_pw[ctt,0,0]
##    array_1[mask] = double_pw[ctt,0,1]
##    array_10[mask] = double_pw[ctt,0,2]
##    array_1m[mask] = double_pw[ctt,0,3]
#    
#    array_sp1[mask] = p_1_single_pw[ctt,0,0]
#    array_sp2[mask] = p_15_single_pw[ctt,0,0]
#    array_sp3[mask] = p_19_single_pw[ctt,0,0]
#    array_sp4[mask] = p_129_single_pw[ctt,0,0]
#    array_dp1[mask] = double_pw[ctt,0,0]
#
#
#new_image_all=np.zeros((620,620,5))
#
#new_image_all[:,:,0]=array_sp1
#new_image_all[:,:,1]=array_sp2
#new_image_all[:,:,2]=array_sp3
#new_image_all[:,:,3]=array_sp4
#new_image_all[:,:,4]=array_dp1


title_all_dust=['dust-0.1um','dust-1um','dust-10um','dust-1mm']
title_all_dust_1=['dust0.1um','dust1um','dust10um','dust1mm']
title_all_pw=['p=1','p=1.5','p=1.9','p=1.29','double p-l']
title_all_band=['4um','5um','8um','24um','70um','160um','870um']

plt.figure(1)
#for ctt in range(7):
for ctt in range(4):
    plt.plot(distance,data_all[:,3,ctt,1],label=str(title_all_dust[ctt]))
#    plt.plot(distance,double_pw[:,ctt,0]/np.median(double_pw[:,ctt,0]),label=str(title_all_band[ctt]))
#    plt.plot(distance,double_pw[:,0,ctt]/np.median(double_pw[:,0,ctt]),label=str(title_all[ctt]))


plt.yscale('log')
plt.xlabel(r'$\rm distance\, (AU)$')
plt.ylabel(r'$\rm L_{\nu}\, [Jy\, cm^{2}]$')
plt.legend()
plt.show()


#np.save('data_all.npy',data_all)

#image_all=generate_image(data_all)
#np.save('image_all.npy',image_all)
image_all=np.load('image_all.npy')


#for ctt in range(7):
#    for ctt2 in range(4):
#        
#        fig, axes = plt.subplots(nrows=2, ncols=3)
#        plt.suptitle(str(title_all_band[ctt])+'  '+str(title_all_dust[ctt2]))
#        
#        for ax, ctt_1 in zip(axes.flat[:-1], range(5)):
#            ax.set_title(str(title_all_pw[ctt_1]))
#        #    ax.imshow(np.log10(np.log10(new_image_all[:,:,ctt]+1)),extent=[-310,310,-310,310])
#            ax.imshow(np.log10(image_all[:,:,ctt,ctt2,ctt_1]),extent=[-310,310,-310,310])
#        
#        ax=axes.flat[-1]
#        for ctt1 in range(5):
#            ax.plot(distance,data_all[:,ctt,ctt2,ctt1],label=str(title_all_pw[ctt1]))
#        ax.set_yscale('log')
#        ax.legend()
#        fig.tight_layout()
#        
#        plt.show()
#        plt.savefig('image_'+str(title_all_band[ctt])+'_'+str(title_all_dust[ctt2])+'.pdf',bbox_inches='tight')
#        plt.clf()
#        plt.close()







"""
xlim=np.array([1,1000])
#ylim=np.array([1e-8,1e15])
#ylim=np.array([1e30,1e50])
ylim=np.array([1e37,10**41.5])
###plot jobs
fig = plt.figure(2)
fig.clf()
ax1 = fig.add_subplot(111)

ax1.plot(lambda_1,dust_spec_01_1,label="dust 0.1um 10AU",color="red",linewidth=2)
ax1.plot(lambda_1,dust_spec_1_1,label="dust 1um 10AU",color="green",linewidth=2)
ax1.plot(lambda_1,dust_spec_10_1,label="dust 10um 10AU",color="blue",linewidth=2)
ax1.plot(lambda_1,dust_spec_1m_1,label="dust 1mm 10AU",color="yellow",linewidth=2)
#
#ax1.plot(lambda_1,dust_spec_01_2,label="dust 0.1um 130AU",color="cyan",linewidth=2)
#ax1.plot(lambda_1,dust_spec_1_2,label="dust 1um 130AU",color="pink",linewidth=2)
#ax1.plot(lambda_1,dust_spec_10_2,label="dust 10um 130AU",color="brown",linewidth=2)
#ax1.plot(lambda_1,dust_spec_1m_2,label="dust 1mm 130AU",color="orange",linewidth=2)

ax1.plot([np.average(fliter_4[:,0]/1e4)],[convole_band_model(dust_spec_01_1,fliter_4)],'o')
ax1.plot([np.average(fliter_5[:,0]/1e4)],[convole_band_model(dust_spec_01_1,fliter_5)],'o')
ax1.plot([np.average(fliter_8[:,0]/1e4)],[convole_band_model(dust_spec_01_1,fliter_8)],'o')
ax1.plot([np.average(fliter_24[:,0]/1e4)],[convole_band_model(dust_spec_01_1,fliter_24)],'o')
ax1.plot([np.average(fliter_70[:,0]/1e4)],[convole_band_model(dust_spec_01_1,fliter_70)],'o')
ax1.plot([np.average(fliter_160[:,0]/1e4)],[convole_band_model(dust_spec_01_1,fliter_160)],'o')
ax1.plot([np.average(fliter_870[:,0]/1e4)],[convole_band_model(dust_spec_01_1,fliter_870)],'o')


ax1.set_xlim(xlim)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.yaxis.get_data_interval()


ax2 = ax1.twiny()
ax2.set_xscale('log')
ax2.set_xlabel(r"$\nu\, [Hz]$",fontsize=15)
ax2.set_xlim(c/xlim*1e4)
#ax2.set_ylim(ylim)

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
#plt.savefig('part_4_fiting_01_sed.pdf',bbox_inches='tight')
plt.show()

"""


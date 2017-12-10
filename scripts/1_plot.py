#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 15:32:49 2017

@author: xuduo
"""


import numpy as np
import matplotlib.pyplot as plt
import astropy
import astropy.constants




fliter_4=np.loadtxt('../filter/Spitzer-IRAC.I2.dat')
fliter_5=np.loadtxt('../filter/Spitzer-IRAC.I3.dat')
fliter_8=np.loadtxt('../filter/Spitzer-IRAC.I4.dat')
fliter_24=np.loadtxt('../filter/Spitzer-MIPS.24mu.dat')
fliter_70=np.loadtxt('../filter/Spitzer-MIPS.70mu.dat')
fliter_160=np.loadtxt('../filter/Spitzer-MIPS.160mu.dat')
fliter_870=np.zeros((10,2))
fliter_870[:,0]=np.linspace(862.8,883,10)*1e4
fliter_870[1:-1,1]=1
title_all_band=['4um','5um','8um','24um','70um','160um','870um']


plt.plot(fliter_4[:,0]/1e4,fliter_4[:,1],label=title_all_band[0])
plt.plot(fliter_5[:,0]/1e4,fliter_5[:,1],label=title_all_band[1])
plt.plot(fliter_8[:,0]/1e4,fliter_8[:,1],label=title_all_band[2])
plt.plot(fliter_24[:,0]/1e4,fliter_24[:,1],label=title_all_band[3])
plt.plot(fliter_70[:,0]/1e4,fliter_70[:,1],label=title_all_band[4])
plt.plot(fliter_160[:,0]/1e4,fliter_160[:,1],label=title_all_band[5])
plt.plot(fliter_870[:,0]/1e4,fliter_870[:,1],label=title_all_band[6])
plt.legend()
plt.xlabel(r"$\lambda\, [\mu m]$",fontsize=15)
plt.ylabel(r'$\rm filter\, response$',fontsize=15)
plt.xscale('log')
plt.savefig('filter.pdf',bbox_inches='tight')

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

temp_db=np.load('temp_double_pw.npy')
temp_1=np.load('temp_p_1_single_pw.npy')
temp_15=np.load('temp_p_1.5_single_pw.npy')
temp_19=np.load('temp_p_1.9_single_pw.npy')
temp_129=np.load('temp_p_1.29_single_pw.npy')

temp_array=np.zeros((200,4,5))
temp_array[:,:,0]=temp_1
temp_array[:,:,1]=temp_15
temp_array[:,:,2]=temp_19
temp_array[:,:,3]=temp_129
temp_array[:,:,4]=temp_db



title_all_dust=['dust-0.1um','dust-1um','dust-10um','dust-1mm']
title_all_dust_1=['dust0.1um','dust1um','dust10um','dust1mm']
title_all_pw=['p=1','p=1.5','p=1.9','p=1.29','double p-l']
title_all_band=['4um','5um','8um','24um','70um','160um','870um']

distance=distance/AU

#plt.clf()
##for ctt1 in range(5):
#for ctt2 in range(4) :
#    plt.plot(distance,temp_array[:,ctt2,0],label=str(title_all_dust_1[ctt2]))
#
#
#plt.xlabel(r'$\rm distance\, (AU)$')
#plt.ylabel(r'$\rm T\, (K)$')
#
#plt.legend()
#plt.show()
#plt.savefig('surface_temp.pdf',bbox_inches='tight')




data_all=np.load('data_all.npy')

plt.clf()
plt.figure(figsize=(10,6))
for ctt1 in range(5):
    for ctt2 in range(4) :
        plt.plot(distance,data_all[:,6,ctt2,ctt1],label=str(title_all_dust_1[ctt2])+' '+str(title_all_pw[ctt1]))


plt.xlabel(r'$\rm distance\, (AU)$')
plt.ylabel(r'$\rm L_{\nu}\, [Jy\, cm^{2}]$')
plt.yscale('log')
plt.title('870 um')
plt.legend()
plt.show()
plt.savefig('870.pdf',bbox_inches='tight')

"""
plt.clf()


for ctt1 in range(7):
    plt.plot(distance,data_all[:,ctt1,0,0],label=str(title_all_band[ctt1]))
plt.title(str(title_all_dust_1[0])+' '+str(title_all_pw[0]))
plt.yscale('log')
plt.xlabel(r'$\rm distance\, (AU)$')
plt.ylabel(r'$\rm L_{\nu}\, [Jy\, cm^{2}]$')
plt.legend()
plt.show()
plt.savefig('distance_band_4um.pdf',bbox_inches='tight')



plt.clf()


data_all=np.load('data_all.npy')


for ctt1 in range(4):
    plt.plot(distance,data_all[:,0,ctt1,0],label=str(title_all_dust_1[ctt1]))
plt.title(str(title_all_band[0])+' '+str(title_all_pw[0]))
plt.yscale('log')
plt.xlabel(r'$\rm distance\, (AU)$')
plt.ylabel(r'$\rm L_{\nu}\, [Jy\, cm^{2}]$')
plt.legend()
plt.show()
plt.savefig('distance_dust_4um.pdf',bbox_inches='tight')

plt.clf()

for ctt1 in range(5):
    plt.plot(distance,data_all[:,6,0,ctt1],label=str(title_all_pw[ctt1]))
plt.title(str(title_all_band[6])+' '+str(title_all_dust_1[0]))
plt.yscale('log')
plt.xlabel(r'$\rm distance\, (AU)$')
plt.ylabel(r'$\rm L_{\nu}\, [Jy\, cm^{2}]$')
plt.legend()
plt.show()
plt.savefig('distance_surfden_870um.pdf',bbox_inches='tight')


"""



image_all=np.load('image_all.npy')
title_all_dust=['dust-0.1um','dust-1um','dust-10um','dust-1mm']
title_all_dust_1=['dust0.1um','dust1um','dust10um','dust1mm']
title_all_pw=['p=1','p=1.5','p=1.9','p=1.29','double p-l']
title_all_band=['4um','5um','8um','24um','70um','160um','870um']



#for ctt in range(7):
#    for ctt2 in range(4):

ctt=3
ctt2 =0
       
fig, axes = plt.subplots(nrows=2, ncols=4,figsize=(20,7))
plt.suptitle(str(title_all_pw[ctt])+'  '+str(title_all_dust[ctt2]))

for ax, ctt_1 in zip(axes.flat[:-1], range(7)):
    ax.set_title(str(title_all_band[ctt_1]))
#    ax.imshow(np.log10(np.log10(new_image_all[:,:,ctt]+1)),extent=[-310,310,-310,310])
    ax.imshow(np.log10(image_all[:,:,ctt_1,ctt2,ctt]),extent=[-310,310,-310,310])

ax=axes.flat[-1]
for ctt1 in range(7):
    ax.plot(distance,data_all[:,ctt1,ctt2,ctt],label=str(title_all_band[ctt1]))
ax.set_yscale('log')
ax.legend()
fig.tight_layout()

plt.show()
plt.savefig('image_'+str(title_all_pw[ctt])+'_'+str(title_all_dust[ctt2])+'.pdf',bbox_inches='tight')
plt.clf()
plt.close()


ctt=4
ctt2 =6
       
#fig, axes = plt.subplots(nrows=2, ncols=4,figsize=(30,10))
fig, axes = plt.subplots(nrows=2, ncols=3)
plt.suptitle(str(title_all_pw[ctt])+'  '+str(title_all_band[ctt2]))

for ax, ctt_1 in zip(axes.flat[:-2], range(4)):
    ax.set_title(str(title_all_dust_1[ctt_1]))
#    ax.imshow(np.log10(np.log10(new_image_all[:,:,ctt]+1)),extent=[-310,310,-310,310])
    ax.imshow(np.log10(image_all[:,:,ctt2,ctt_1,ctt]),extent=[-310,310,-310,310])

ax=axes.flat[-2]
ax = plt.subplot2grid((2, 3), (1, 1), colspan=2)

for ctt1 in range(4):
    ax.plot(distance,data_all[:,ctt2,ctt1,ctt],label=str(title_all_dust_1[ctt1]))
ax.set_yscale('log')
ax.legend()
fig.tight_layout()

plt.show()
plt.savefig('image_'+str(title_all_pw[ctt])+'_'+str(title_all_band[ctt2])+'.pdf',bbox_inches='tight')
plt.clf()
plt.close()


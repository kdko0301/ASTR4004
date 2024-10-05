# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:41:36 2024

@author: kdko0
"""

#imports
import illustris_python as il
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py


def Plot_2D(num):
    ### Create of 2D plot for the FOF halo system with the given halo index
    
    #Import data
    basePath = './TNG100' #set basepath

    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'Group_R_Crit200', 'GroupNsubs', 'GroupCM'] #For FOF halos
    sub_fields = [ 'SubhaloCM', 'SubhaloMass', 'SubhaloHalfmassRadType'] #For subfind subhalos

    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)


    FOF_CM= halos['GroupCM'][num]
    FOF_R200 = halos['Group_R_Crit200'][num]
    sub_count = halos['GroupNsubs'][num]
    pri_num = halos['GroupFirstSub'][num]
    pri_CM = subhalos['SubhaloCM'][pri_num]
    pri_r_half = subhalos['SubhaloHalfmassRadType'][pri_num][4]
    #pri_R200 = pri_r_half/0.015 # From the half stellar mass radius - R200 relation
    pri_R200 = R_2t_vdis(pri_num)
    pri_R200_mass = Mass_R200(pri_num)
    pri_R_2t = R_2t_mass(pri_num)
    sub_x_coords = []
    sub_y_coords = []
    lim = max(FOF_R200, pri_R200, pri_R200_mass, pri_R_2t)
    diff = max(abs(FOF_CM[0]-pri_CM[0]), abs( FOF_CM[1]-pri_CM[1]))
    theta = np.linspace( 0 , 2 * np.pi , 150 )
    for i in range(pri_num+1, pri_num+sub_count):
        sub_x_coords.append(subhalos['SubhaloCM'][i][0]-FOF_CM[0])
        sub_y_coords.append(subhalos['SubhaloCM'][i][1]-FOF_CM[1])
    print ('Halo index:', num,'Subhalo count:', sub_count)
    print('FOF R200:', FOF_R200, 'Primary subhalo R200:', pri_R200)

    plt.plot(sub_x_coords, sub_y_coords, 'c.', alpha = 0.5)
    plt.plot (0,0, 'k*')
    plt.plot(pri_CM[0]-FOF_CM[0], pri_CM[1]-FOF_CM[1], 'r*')
    plt.plot(FOF_R200 * np.cos( theta ), FOF_R200 * np.sin( theta ), 'b-')
    plt.plot(pri_CM[0]-FOF_CM[0] + pri_R200 * np.cos( theta ),pri_CM[1]-FOF_CM[1] + pri_R200 * np.sin( theta ), 'm-')
    plt.plot(pri_CM[0]-FOF_CM[0] + pri_R200_mass * np.cos( theta ),pri_CM[1]-FOF_CM[1] + pri_R200_mass * np.sin( theta ), 'g-')
    plt.plot(pri_CM[0]-FOF_CM[0] + pri_R_2t * np.cos( theta ),pri_CM[1]-FOF_CM[1] + pri_R_2t * np.sin( theta ), 'y-')
    plt.xlim(-lim -diff, lim +diff)
    plt.ylim(-lim-diff, lim+diff)
    plt.title('halo number: {num}, subhalo count: {sub_count}'.format(num = num, sub_count = sub_count))
    plt.xlabel('Relative x-coordinates (ckpc/h)')
    plt.ylabel('Relative y-coordinates (ckpc/h)')
    plt.legend(['subhalos CM', 'FOF CM', 'Primary CM', 'FOF R200', 'Size R200', 'Mass R200', 'R_2t'], loc='upper right', fontsize = 'x-small')
    
    return None
def Sort(low_lim, high_lim):
    ### Give index of halos with the subhalo count within the given range 
    basePath = './TNG100' #set basepath

    #Set fields for data import
    FOF_fields = [ 'GroupNsubs', 'GroupFirstSub' ] #For FOF halos

    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)

    ### Return index of FOF halos with subhalo count in the given range
    num_list = []
    i = 0
    for count in halos['GroupNsubs']:
        if low_lim <= count <= high_lim:
            num_list.append(i)
        i+=1
    print('Counts:', len(num_list))
    return num_list
    
def Isolation(max_dist, index):
    ### Return the distance to the cloesest FOFhalo within the max distance and max index
    basePath = './TNG100' #set basepath
    #Set fields for data import
    FOF_fields = [ 'GroupCM'] #For FOF halos
    #data import from TNG100 snapshot 99 (z=0) groupcat
    FOF_CM = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)




    iso = max_dist
    CM = FOF_CM[index]
    for CM2 in FOF_CM:
        distance = np.sqrt((CM2[0]-CM[0])**2+(CM2[1]-CM[1])**2+(CM2[2]-CM[2])**2)
        if distance < 10:
            pass
        elif distance < iso:
            iso = distance
        
    return iso

def Mass_R200(index):
    ### Calculate R200 based on assumption that total mass equal M200
    
    #Import data
    basePath = './TNG100' #set basepath
    sub_fields = [ 'SubhaloCM', 'SubhaloMass']
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)
    sub_mass = subhalos['SubhaloMass'][index]
    
    p_c = 125.2e-10 #critical density in units of 10^10 solar mass per kpc^3
    R200 = (sub_mass/(4/3*np.pi*200*p_c))**(1/3)
    return R200

def Size_R200(index):
    ### Calculate R200 based on the half stellar mass radius
    #Import data
    basePath = './TNG100' #set basepath

    #Set fields for data import

    sub_fields = [ 'SubhaloCM', 'SubhaloHalfmassRadType'] #For subfind subhalos

    #data import from TNG100 snapshot 99 (z=0) groupcat

    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)


    sub_r_half = subhalos['SubhaloHalfmassRadType'][index][4]
    R200 = sub_r_half/0.015 # From the half stellar mass radius - R200 relation
    return R200

def R_2t_mass(index):
    ### Calculate R_2t based on the subhalo mass
    
    #Import data
    basePath = './TNG100' #set basepath
    sub_fields = [  'SubhaloMass']
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)
    sub_mass = subhalos[index]
    

    R_2t = 0.215*(sub_mass/100)**(1/3)*1000 #R_2t in units of kpc
    return R_2t

def R_2t_vdis(index):
    ### Calculate R_2t based on the velocity dispersion of galaxy elements
    
    #Import data
    basePath = './TNG100' #set basepath
    sub_fields = [ 'SubhaloVelDisp']
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)
    sub_vdis = subhalos[index]   
    R_2t = sub_vdis/368*1000
    return R_2t

def Halo_Condition(index):
    ### Consider FOF halo conditions and verify whether it is insiginificant or not
    #Import data
    basePath = './TNG100' #set basepath

    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'GroupNsubs', 'GroupLenType', 'GroupMassType'] 
    sub_fields = [ 'SubhaloFlag', 'SubhaloHalfmassRadType', 'SubhaloStellarPhotometrics']
    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)
    
    # Primary subhalo flag 
    pri_num = halos['GroupFirstSub'][index]
    flag = subhalos['SubhaloFlag'][pri_num]
    star_count = halos['GroupLenType'][index][4]
    star_mass = halos['GroupMassType'][index][4]
    M_r = subhalos['SubhaloStellarPhotometrics'][pri_num][5]
    if flag < 1:
        Cond = 0
    # Star existense
    elif M_r >-18:
        Cond = 0
    
    elif star_count < 100:
        Cond = 0
    # Star Mass
    
    elif star_mass < 1e-7:
        Cond = 0
    
    else:
        Cond = 1
    
    return Cond

def Halo_list():
    ### Get index list of halos that satisfy the conditions
    
    basePath = './TNG100' #set basepath

    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'GroupNsubs', 'GroupLenType', 'GroupMassType'] 
    sub_fields = [ 'SubhaloFlag', 'SubhaloMass','SubhaloMassType', 'SubhaloStellarPhotometrics']
    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)
    halocount = halos['count']
    halo_list = []
    
    pri_num = halos['GroupFirstSub']
    flag_list = subhalos['SubhaloFlag']
    Len = halos['GroupLenType']
    Mass = halos['GroupMassType']
    Mag = subhalos['SubhaloStellarPhotometrics']
    sub_mass = subhalos['SubhaloMass']
    counter = 0

    for i in range(0, halocount-1):
        pri = pri_num[i]
        flag = flag_list[pri]

        M_r = Mag[pri][5]
        sub_count = halos['GroupNsubs'][i]
        pri_mass = subhalos['SubhaloMass'][pri]
        
        if flag > 0:
        # Star existense
            if M_r <-18 and pri_mass <1.5e3 and pri_mass > 1.5e1:
                if sub_count>4:
                    
                    if Mag[pri+1][5] <-10 and Mag[pri+2][5] <-10:

                        halo_list.append(i)
        counter += 1
        if counter == 1000000:
            print('Listing Progress:', i/halocount*100, '%' )
            counter = 0
    
    return halo_list

def Halo_list_rv():
    ### Get index list of halos that satisfy the conditions
    
    basePath = './TNG100' #set basepath

    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'GroupNsubs', 'GroupLenType', 'GroupMassType'] 
    sub_fields = [ 'SubhaloFlag', 'SubhaloMass','SubhaloMassType', 'SubhaloStellarPhotometrics']
    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)
    halocount = halos['count']
    halo_list = []
    
    pri_num = halos['GroupFirstSub']
    flag_list = subhalos['SubhaloFlag']
    Len = halos['GroupLenType']
    Mass = halos['GroupMassType']
    Mag = subhalos['SubhaloStellarPhotometrics']
    sub_mass = subhalos['SubhaloMass']
    counter = 0

    for i in range(0, halocount-1):
        pri = pri_num[i]
        flag = flag_list[pri]
        star_mass = Mass[i][4]
        star_count = Len[i][4]
        M_r = Mag[pri][5]
        sub_count = halos['GroupNsubs'][i]
        pri_mass = subhalos['SubhaloMass'][pri]
        
        if flag > 0:
        # Star existense
            if M_r <-18 and pri_mass <1.5e3 and pri_mass > 1.5e1:
                if sub_count>4:
                    sec_mass = sub_mass[pri+1]
                    tri_mass = sub_mass[pri+2]
                    if sec_mass< 0.2 and tri_mass < 0.2:

                        halo_list.append(i)
        counter += 1
        if counter == 1000000:
            print('Listing Progress:', i/halocount*100, '%' )
            counter = 0
    
    return halo_list
def Halo_list_nl():
    ### Get index list of halos that satisfy the conditions
    
    basePath = './TNG100' #set basepath

    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'GroupNsubs', 'GroupLenType', 'GroupMassType'] 
    sub_fields = [ 'SubhaloFlag', 'SubhaloMass', 'SubhaloMassType', 'SubhaloStellarPhotometrics']
    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)
    halocount = halos['count']
    halo_list = []
    
    pri_num = halos['GroupFirstSub']
    flag_list = subhalos['SubhaloFlag']
    Len = halos['GroupLenType']
    Mass = halos['GroupMassType']
    Mag = subhalos['SubhaloStellarPhotometrics']
    sub_mass = subhalos['SubhaloMassType']
    
    counter = 0

    for i in range(0, halocount-1):
        pri = pri_num[i]
        flag = flag_list[pri]
        star_mass = Mass[i][4]
        star_count = Len[i][4]
        M_r = Mag[pri][5]
        sub_count = halos['GroupNsubs'][i]
        pri_mass = subhalos['SubhaloMass'][pri]
        if flag > 0:
        # Star existense
            if M_r <-18 and pri_mass <1.5e3 and pri_mass > 1.5e1:
                if sub_count>4:
                    halo_list.append(i)
        counter += 1
        if counter == 1000000:
            print('Listing Progress:', i/halocount*100, '%' )
            counter = 0
    
    return halo_list
def Save_halo_list(filename):
    halo_list = Halo_list()
    dict = {'Halo_ind': halo_list}

    df = pd.DataFrame(dict)
    df.to_csv(filename)
    return None
    
def Evaluate_nolim(index):
    #Evaluate the galaxy group paramters with no r2t, halo condition limit
    
    #Import data
    basePath = './TNG100' #set basepath

    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'GroupNsubs', 'GroupPos'] #For FOF halos
    sub_fields = [ 'SubhaloFlag' ,'SubhaloCM', 'SubhaloMass', 'SubhaloVelDisp', 'SubhaloStellarPhotometrics'] #For subfind subhalos
    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)    
    
    # Get required variables
    sub_count = halos['GroupNsubs'][index]
    pri_num = halos['GroupFirstSub'][index]
    pri_CM = subhalos['SubhaloCM'][pri_num]
    pri_R2t = R_2t_vdis(pri_num)
    pri_R200 = Mass_R200(pri_num)
    
    # Obtain subhalo Coords and delta_M
    sub_x_coords = []
    sub_y_coords = []
    pri_mag = subhalos['SubhaloStellarPhotometrics'][pri_num]
    pri_flux = 10**(-pri_mag[4]/2.5) + 10**(-pri_mag[5]/2.5) + 10**(-pri_mag[6]/2.5) +10**(-pri_mag[7]/2.5) # in units of 3631Jy
    flux_diff_list = []
    for i in range(pri_num+1, pri_num+sub_count):
        sub_x_coords.append(subhalos['SubhaloCM'][i][0]-pri_CM[0])
        sub_y_coords.append(subhalos['SubhaloCM'][i][1]-pri_CM[1])
        sub_mag = subhalos['SubhaloStellarPhotometrics'][i]
        sub_flux = 10**(-sub_mag[4]/2.5) + 10**(-sub_mag[5]/2.5) + 10**(-sub_mag[6]/2.5) +10**(-sub_mag[7]/2.5) # in units of 3631Jy
        flux_diff = pri_flux - sub_flux
        flux_diff_list.append(flux_diff)
    
    sort_diff_list = sorted(flux_diff_list)
    sec_ind = flux_diff_list.index(sort_diff_list[0])
    tri_ind = flux_diff_list.index(sort_diff_list[1])
    sec_num = pri_num+1+sec_ind
    tri_num = pri_num+1+tri_ind
    sec_mag = subhalos['SubhaloStellarPhotometrics'][sec_num]
    tri_mag = subhalos['SubhaloStellarPhotometrics'][tri_num]
    mag_diff_sec = np.abs(pri_mag-sec_mag)
    mag_diff_tri = np.abs(pri_mag-tri_mag)
    
    # Get the distance to the closest halo
    FOF_pos = halos['GroupPos'][index]
    dist_list = []

    for pos in halos['GroupPos']:
        dist = np.sqrt((pos[0]-FOF_pos[0])**2+(pos[1]-FOF_pos[1])**2+(pos[2]-FOF_pos[2])**2)
        dist_list.append(dist)
    sort_dist_list = sorted(dist_list)
    min_dist = sort_dist_list[1]
    if min_dist < pri_R2t:
        iso = 'False'
    else:
        iso = 'True'
    
    # Create 2D- projected  plot
    sec_CM = subhalos['SubhaloCM'][sec_num]
    tri_CM = subhalos['SubhaloCM'][tri_num]
    theta = np.linspace( 0 , 2 * np.pi , 150 )
    max_lim = max(max(sub_x_coords), max(sub_y_coords), pri_R2t)
    min_lim = min(min(sub_x_coords), min(sub_y_coords), -pri_R2t)
    plt.figure(figsize = (8, 8))
    plt.xlim(min_lim, max_lim)
    plt.ylim(min_lim, max_lim)
    plt.plot(sub_x_coords, sub_y_coords, 'k.', alpha = 0.5, ms = 5)
    plt.plot(0, 0, 'bo', ms = 10)
    plt.plot(sec_CM[0]-pri_CM[0], sec_CM[1]-pri_CM[1], 'go', ms = 8)
    plt.plot(tri_CM[0]-pri_CM[0], tri_CM[1]-pri_CM[1], 'ro', ms = 8)
    plt.plot(pri_R2t * np.cos( theta ), pri_R2t * np.sin( theta ), 'b--', ms = 5)
    plt.title('halo number: {index}, subhalo count: {sub_count}'.format(index = index, sub_count = sub_count))
    plt.xlabel('Relative x-coordinates (ckpc/h)')
    plt.ylabel('Relative y-coordinates (ckpc/h)')
    plt.legend(['subhalos CM', 'Primary CM', 'Secondary CM', 'Tertiary CM', 'R_2t'], loc='upper right', fontsize = 'x-small')
    
    # Print galaxy parameters
    print('FOF index:', index, 'Subhalo count', sub_count)
    print('R_2t(ckpc/h):', pri_R2t, 'R200(ckpc/h):', pri_R200, 'Min dist(ckpc/h):', min_dist, 'Isolation:', iso )
    print('Primary index:', pri_num, 'Secondary index:', sec_num, 'Tertiary index:', tri_num)
    print('Delta M (secondary):', 'g:', mag_diff_sec[4], 'r:', mag_diff_sec[5], 'i:', mag_diff_sec[6], 'z:', mag_diff_sec[7] )
    print('Delta M (Tertiary):', 'g:', mag_diff_tri[4], 'r:', mag_diff_tri[5], 'i:', mag_diff_tri[6], 'z:', mag_diff_tri[7] )
    
    return mag_diff_sec, mag_diff_tri

def Evaluate(index):
    #Evaluate the galaxy group paramters with r2t, halo condition limit
    
    #Import data
    basePath = './TNG100' #set basepath

    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'GroupNsubs', 'GroupPos'] #For FOF halos
    sub_fields = [ 'SubhaloFlag' ,'SubhaloCM', 'SubhaloMass', 'SubhaloVelDisp', 'SubhaloStellarPhotometrics'] #For subfind subhalos
    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)    
    
    #Test for condtions
    if Halo_Condition(index) < 1:
        print('Insignificant halo')
        return None
    
    # Get required variables
    sub_count = halos['GroupNsubs'][index]
    pri_num = halos['GroupFirstSub'][index]
    pri_CM = subhalos['SubhaloCM'][pri_num]
    pri_R2t = subhalos['SubhaloVelDisp'][pri_num]/368*1000
    pri_R200 = Mass_R200(pri_num)
    
    # Obtain subhalo Coords and delta_M
    sub_x_coords = []
    sub_y_coords = []
    pri_mag = subhalos['SubhaloStellarPhotometrics'][pri_num]
    pri_mass = subhalos['SubhaloMass'][pri_num]
    pri_flux = 10**(-pri_mag[4]/2.5) + 10**(-pri_mag[5]/2.5) + 10**(-pri_mag[6]/2.5) +10**(-pri_mag[7]/2.5) # in units of 3631Jy
    flux_diff_list = []
    sat_count = 0
    for i in range(pri_num+1, pri_num+sub_count):
        sub_mag = subhalos['SubhaloStellarPhotometrics'][i]
        if sub_mag[4] < -10:
            sub_x_coord = subhalos['SubhaloCM'][i][0]-pri_CM[0]
            sub_y_coord = subhalos['SubhaloCM'][i][1]-pri_CM[1]
            sub_z_coord = subhalos['SubhaloCM'][i][2]-pri_CM[2]
            sub_x_coords.append(sub_x_coord)
            sub_y_coords.append(sub_y_coord)
            dist = np.sqrt(sub_x_coord**2+sub_y_coord**2+sub_z_coord**2)
            if dist > pri_R2t:
                flux_diff_list.append(1e99) #Arbitary large value (non realistic)
            else:
                sub_mag = subhalos['SubhaloStellarPhotometrics'][i]
                sub_flux = 10**(-sub_mag[4]/2.5) + 10**(-sub_mag[5]/2.5) + 10**(-sub_mag[6]/2.5) +10**(-sub_mag[7]/2.5) # in units of 3631Jy
                flux_diff = pri_flux - sub_flux
                flux_diff_list.append(flux_diff)
                sat_count+=1
    sort_diff_list = sorted(flux_diff_list)

    sec_ind = flux_diff_list.index(sort_diff_list[0])
    tri_ind = flux_diff_list.index(sort_diff_list[1])
    sec_num = pri_num+1+sec_ind
    tri_num = pri_num+1+tri_ind
    sec_mag = subhalos['SubhaloStellarPhotometrics'][sec_num]
    tri_mag = subhalos['SubhaloStellarPhotometrics'][tri_num]
    sec_mass =  subhalos['SubhaloMass'][sec_num]
    tri_mass =  subhalos['SubhaloMass'][tri_num]
    mag_diff_sec = np.abs(pri_mag-sec_mag)
    mag_diff_tri = np.abs(pri_mag-tri_mag)
    
    
    # Create 2D- projected  plot
    sec_CM = subhalos['SubhaloCM'][sec_num]
    tri_CM = subhalos['SubhaloCM'][tri_num]
    theta = np.linspace( 0 , 2 * np.pi , 150 )
    max_lim = 1.5 *pri_R2t
    min_lim = -1.5*pri_R2t
    plt.figure(figsize = (8, 8))
    plt.xlim(min_lim, max_lim)
    plt.ylim(min_lim, max_lim)
    plt.plot(sub_x_coords, sub_y_coords, 'k.', alpha = 0.5, ms = 5)
    plt.plot(0, 0, 'bo', ms = 10)
    plt.plot(sec_CM[0]-pri_CM[0], sec_CM[1]-pri_CM[1], 'go', ms = 8)
    plt.plot(tri_CM[0]-pri_CM[0], tri_CM[1]-pri_CM[1], 'ro', ms = 8)
    plt.plot(pri_R2t * np.cos( theta ), pri_R2t * np.sin( theta ), 'b--', ms = 5)
    plt.title('halo number: {index}, subhalo count: {sub_count}, satelite count: {sat_count}'.format(index = index, sub_count = sub_count, sat_count = sat_count))
    plt.xlabel('Relative x-coordinates (ckpc/h)')
    plt.ylabel('Relative y-coordinates (ckpc/h)')
    plt.legend(['subhalos CM', 'Primary CM', 'Secondary CM', 'Tertiary CM', 'R_2t'], loc='upper right', fontsize = 'x-small')
    
    # Print galaxy parameters
    print('FOF index:', index, 'Subhalo count', sub_count, 'satelite count:', sat_count)
    print('R_2t(ckpc/h):', pri_R2t, 'R200(ckpc/h):', pri_R200)
    print('Primary index:', pri_num, 'Secondary index:', sec_num, 'Tertiary index:', tri_num)
    print('Primary mass:', pri_mass, 'Secondary mass:', sec_mass, 'Tertiary mass:', tri_mass)
    print('Delta M (secondary):', 'g:', mag_diff_sec[4], 'r:', mag_diff_sec[5], 'i:', mag_diff_sec[6], 'z:', mag_diff_sec[7] )
    print('Delta M (Tertiary):', 'g:', mag_diff_tri[4], 'r:', mag_diff_tri[5], 'i:', mag_diff_tri[6], 'z:', mag_diff_tri[7] )
    
    return None

def Full_Evaluate(halo_list):
    #Evaluate the galaxy group paramters for list of halos
    
    #Import data
    basePath = './TNG100' #set basepath

    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'GroupNsubs', 'GroupPos'] #For FOF halos
    sub_fields = [ 'SubhaloFlag' ,'SubhaloCM', 'SubhaloMass', 'SubhaloVelDisp', 'SubhaloStellarPhotometrics'] #For subfind subhalos
    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)  
    
    # For progress report
    max_count = len(halo_list)
    progress = 0
    counter = 0
    # Set empty lists
    iso_halo_list = []
    delta_M_sec = []
    delta_M_tri = []
    sec_dist_list = []
    tri_dist_list = []
    sub_count_list = []
    R2t_list = []
    vdis_list = []
    pri_mag_list = []
    delta_mass_sec = []
    delta_mass_tri = []
    pri_mass_list = []
    # Evaluate and Append
    for index in halo_list:
        # Progress meter
        progress += 1
        counter += 1
        if counter == 1000:
            print('Progress:', progress/max_count*100, '%' )
            counter = 0
        
        # Get base variables
        sub_count = halos['GroupNsubs'][index]
        pri_num = halos['GroupFirstSub'][index]
        pri_CM = subhalos['SubhaloCM'][pri_num]
        pri_vdis = subhalos['SubhaloVelDisp'][pri_num]
        pri_R2t = pri_vdis/368*1000
        # Get the distance to the closest halo
        FOF_pos = halos['GroupPos'][index]
        dist_list = []
        for i in halo_list:
            pos = halos['GroupPos'][i]
            dist = np.sqrt((pos[0]-FOF_pos[0])**2+(pos[1]-FOF_pos[1])**2+(pos[2]-FOF_pos[2])**2)
            dist_list.append(dist)
        sort_dist_list = sorted(dist_list)
        min_dist = sort_dist_list[1]
        if min_dist < pri_R2t:
            pass
        else:       
            # Obtain subhalo Coords and delta_M
            sub_dist_list = []
            flux_diff_list = []
            mass_diff_list = []
            pri_mag = subhalos['SubhaloStellarPhotometrics'][pri_num]
            pri_flux = 10**(-pri_mag[4]/2.5) + 10**(-pri_mag[5]/2.5) + 10**(-pri_mag[6]/2.5) +10**(-pri_mag[7]/2.5) # in units of 3631Jy
            pri_mass = subhalos['SubhaloMass'][pri_num]
            
            for i in range(pri_num+1, pri_num+sub_count):
                sub_x_coord = subhalos['SubhaloCM'][i][0]-pri_CM[0]
                sub_y_coord = subhalos['SubhaloCM'][i][1]-pri_CM[1]
                sub_z_coord = subhalos['SubhaloCM'][i][2]-pri_CM[2]

                sub_dist = np.sqrt(sub_x_coord**2+sub_y_coord**2+sub_z_coord**2)
                sub_dist_list.append(dist)
                if sub_dist > pri_R2t:
                    flux_diff_list.append(1e99) #Arbitary large value (non realistic)
                    mass_diff_list.append(1e99)
                else:
                    sub_mag = subhalos['SubhaloStellarPhotometrics'][i]
                    sub_flux = 10**(-sub_mag[4]/2.5) + 10**(-sub_mag[5]/2.5) + 10**(-sub_mag[6]/2.5) +10**(-sub_mag[7]/2.5) # in units of 3631Jy
                    flux_diff = pri_flux - sub_flux
                    flux_diff_list.append(flux_diff)
                    sub_mass = subhalos['SubhaloMass'][i]
                    mass_diff = np.abs(pri_mass - sub_mass)
                    mass_diff_list.append(mass_diff)
            sort_diff_list = sorted(flux_diff_list)
            sec_ind = flux_diff_list.index(sort_diff_list[0])
            tri_ind = flux_diff_list.index(sort_diff_list[1])
            sec_num = pri_num+1+sec_ind
            tri_num = pri_num+1+tri_ind
            sec_mag = subhalos['SubhaloStellarPhotometrics'][sec_num]
            tri_mag = subhalos['SubhaloStellarPhotometrics'][tri_num]
            mag_diff_sec = np.abs(pri_mag-sec_mag)
            mag_diff_tri = np.abs(pri_mag-tri_mag)
            sec_dist = sub_dist_list[sec_ind]
            tri_dist = sub_dist_list[tri_ind]
            
            sort_mass_list = sorted(mass_diff_list)
            sec_ind_mass = mass_diff_list.index(sort_mass_list[0])
            tri_ind_mass = mass_diff_list.index(sort_mass_list[1])
            sec_num_mass = pri_num+1+sec_ind_mass
            tri_num_mass = pri_num+1+tri_ind_mass
            sec_mass = subhalos['SubhaloMass'][sec_num_mass]
            tri_mass = subhalos['SubhaloMass'][tri_num_mass]
            mass_diff_sec = np.abs(np.log10(pri_mass) - np.log10(sec_mass))
            mass_diff_tri = np.abs(np.log10(pri_mass) - np.log10(tri_mass))
            
            if mag_diff_sec[5]<1e36:
                iso_halo_list.append(index)
                delta_M_sec.append(mag_diff_sec)
                delta_M_tri.append(mag_diff_tri)
                sec_dist_list.append(sec_dist)
                tri_dist_list.append(tri_dist)
                sub_count_list.append(sub_count)
                R2t_list.append(pri_R2t)
                vdis_list.append(pri_vdis)
                pri_mag_list.append(pri_mag)
                delta_mass_sec.append(mass_diff_sec)
                delta_mass_tri.append(mass_diff_tri)
                pri_mass_list.append(pri_mass)
            else:
                pass
            
    return iso_halo_list, delta_M_sec, delta_M_tri, sec_dist_list, tri_dist_list, sub_count_list, R2t_list, vdis_list, pri_mag_list, delta_mass_sec, delta_mass_tri, pri_mass_list
        
def GLF(index, bin_num):
    #Evaluate the galaxy group paramters with r2t, halo condition limit

    #Import data
    basePath = './TNG100' #set basepath

    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'GroupNsubs', 'GroupPos'] #For FOF halos
    sub_fields = [ 'SubhaloFlag' ,'SubhaloCM', 'SubhaloMass','SubhaloMassType', 'SubhaloVelDisp', 'SubhaloStellarPhotometrics'] #For subfind subhalos
    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)            
    sub_count = halos['GroupNsubs'][index]
    pri_num = halos['GroupFirstSub'][index]
    pri_CM = subhalos['SubhaloCM'][pri_num]
    pri_vdis = subhalos['SubhaloVelDisp'][pri_num]
    pri_R2t = pri_vdis/368*1000
    pri_mag = subhalos['SubhaloStellarPhotometrics'][pri_num][4]

   
     # Obtain subhalo Coords and delta_M

    mag_list = []
    mass_list = []

     
    for i in range(pri_num, pri_num+sub_count):
        sub_x_coord = subhalos['SubhaloCM'][i][0]-pri_CM[0]
        sub_y_coord = subhalos['SubhaloCM'][i][1]-pri_CM[1]
        sub_z_coord = subhalos['SubhaloCM'][i][2]-pri_CM[2]

        sub_dist = np.sqrt(sub_x_coord**2+sub_y_coord**2+sub_z_coord**2)
        sub_mag = subhalos['SubhaloStellarPhotometrics'][i][4]
        sub_mass = subhalos['SubhaloMassType'][i][1]
        if sub_dist > pri_R2t:
            pass
        elif sub_mag>-10:
            pass
        elif sub_mass < 1e-10:
            pass
        else:
            
            mag_list.append(sub_mag)
            mass_list.append(sub_mass)
    print( "subhalo count:", len(mag_list))
    
    plt.figure(0)
    fig, axs = plt.subplots(2)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
    fig.suptitle('GLF and DMMF')
    axs[0].hist(mag_list, bins = bin_num)

    axs[0].set_xlabel('$ M_r $')
    axs[0].set_ylabel('Frequency')
    
    axs[1].hist(np.log10(mass_list), bins = bin_num)

    axs[1].set_xlabel('$log_{10}$(DM mass) ($log_{10}(M_{sol}*10^{10}))$')
    axs[1].set_ylabel('Frequency')
    plt.gca().invert_xaxis()
 
    return None
    
def L_gal_sublist():
    basePath = './TNG100' #set basepath
    path = './TNG100/L_Galaxies/LGalaxies_099.hdf5'

    FOF_fields = ['GroupFirstSub', 'GroupNsubs'] #For FOF halos
    sub_fields = [  'SubhaloMass', 'SubhaloVelDisp', 'SubhaloStellarPhotometrics', 'SubhaloGrNr'] #For subfind subhalos

    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)

    TNGMag = subhalos['SubhaloStellarPhotometrics']
    submass = subhalos['SubhaloMass']
    #halo_list = fc.Halo_list()
    # Import data
    L_Gal=  h5py.File(path, 'r')['Galaxy']

    MagDust = np.array(L_Gal['MagDust'])

    TNG_index = np.array(L_Gal['SubhaloIndex_TNG'])
    Type = np.array(L_Gal['Type'])


    index_list = []
    FOF_list = []


    i = 0
    while i < len(TNG_index):
        if TNG_index[i]>=0:
            if Type[i] == 0:
                if MagDust[i][17] <-18:
                    if submass[TNG_index[i]] <1.5e3 and submass[TNG_index[i]] >1.5e1:
                        index_list.append(TNG_index[i])
                        FOF_list.append(subhalos['SubhaloGrNr'][TNG_index[i]])          
                        
        i+=1
    return index_list, FOF_list

def Evaluate_LGal(index):
    #Evaluate the galaxy group paramters with r2t, halo condition limit
    
    #Import data
    basePath = './TNG100' #set basepath
    path = './TNG100/L_Galaxies/LGalaxies_099.hdf5'
    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'GroupNsubs', 'GroupPos'] #For FOF halos
    sub_fields = [ 'SubhaloGrNr', 'SubhaloFlag' ,'SubhaloCM', 'SubhaloMass','SubhaloMassType', 'SubhaloVelDisp', 'SubhaloStellarPhotometrics'] #For subfind subhalos
    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)    
    L_Gal=  h5py.File(path, 'r')['Galaxy']

    MagDust = np.array(L_Gal['MagDust'])
    TNG_index = np.array(L_Gal['SubhaloIndex_TNG'])

    CentDist = np.array(L_Gal['DistanceToCentralGal'])
    StellarMass = np.array(L_Gal['StellarMass'])
    XrayLum = np.array(L_Gal['XrayLum'])
    vel_dis = subhalos['SubhaloVelDisp']
    R2t = vel_dis/368*1000
    #Test for condtions
    
    
    pri = TNG_index[index]
    FOF  = subhalos['SubhaloGrNr'][pri]
    sub_count = halos['GroupNsubs'][FOF]
    pri_mag = MagDust[index]
    pri_mass = subhalos['SubhaloMass'][pri]
    pri_DMmass = subhalos['SubhaloMassType'][pri][1]
    pri_starmass = StellarMass[index]


    pri_R2t = R2t[pri]
    delta_M_list = []
    delta_DM_list = []
    delta_star_list = []
    delta_mass_list =[]
    mass_list =[]
    SM_list =[]
    x_list = []
    y_list =[]
    count = 0

    for sub in range(pri+1, pri+ sub_count):

        if sub in TNG_index:
            L_ind  = np.where(TNG_index == sub)[0][0]

            sub_mag = MagDust[L_ind]
            
            dist  = np.sqrt(CentDist[L_ind][0]**2+CentDist[L_ind][1]**2+CentDist[L_ind][2]**2)
            if dist<= pri_R2t and sub_mag[16]<99:
                count += 1
                sub_mass = subhalos['SubhaloMass'][sub]
                sub_DMmass = subhalos['SubhaloMassType'][sub][1]
                sub_starmass = StellarMass[L_ind]
                mag_diff  = np.abs(pri_mag[16]-sub_mag[16])
                mass_diff = np.abs(np.log10(pri_mass)-np.log10(sub_mass))
                starmass_diff = np.abs(np.log10(pri_starmass)-np.log10(sub_starmass))
                DM_diff = np.abs(np.log10(pri_DMmass)-np.log10(sub_DMmass))
                delta_M_list.append(mag_diff)
                delta_DM_list.append(DM_diff)
                delta_star_list.append(starmass_diff)
                delta_mass_list.append(mass_diff)
                sub_x  = CentDist[L_ind][0]
                sub_y  = CentDist[L_ind][1]
                x_list.append(sub_x)
                y_list.append(sub_y)
                mass_list.append(subhalos['SubhaloMass'][sub])
                SM_list.append(StellarMass[L_ind])
    if count >1:
        delta_M = sorted(delta_M_list)[0]
        delta_DM = sorted(delta_DM_list)[0]
        delta_star = sorted(delta_star_list)[0]
        delta_mass = sorted(delta_mass_list)[0]
        sec = delta_M_list.index(delta_M)
        sec_x = x_list[sec]
        sec_y = y_list[sec]
        sec_mass = mass_list[sec]
        sec_SM = SM_list[sec]
    else:
        print('isolated galaxy')



    
    # Create 2D- projected  plot

    theta = np.linspace( 0 , 2 * np.pi , 150 )
    max_lim = 1.5 *pri_R2t
    min_lim = -1.5*pri_R2t
    plt.figure(figsize = (8, 8))
    plt.xlim(min_lim, max_lim)
    plt.ylim(min_lim, max_lim)
    plt.plot(x_list, y_list, 'k.', alpha = 0.5, ms = 5)
    plt.plot(0, 0, 'bo', ms = 10)
    plt.plot(sec_x, sec_y, 'go', ms = 8)
   
    plt.plot(pri_R2t * np.cos( theta ), pri_R2t * np.sin( theta ), 'b--', ms = 5)
    plt.title('Galaxy number: {index}, subhalo count: {sub_count}, galaxy count: {count}'.format(index = index, sub_count = sub_count, count = count))
    plt.xlabel('Relative x-coordinates (ckpc/h)')
    plt.ylabel('Relative y-coordinates (ckpc/h)')
    plt.legend(['subhalos CM', 'Primary CM', 'Secondary CM',  'R_2t'], loc='upper right', fontsize = 'x-small')
    
    # Print galaxy parameters
    print('FOF index:', FOF, 'Subhalo count', sub_count, 'galaxy count:', count)
    print('R_2t(ckpc/h):', pri_R2t)
    print('Primary mass:', pri_mass, 'Secondary mass:', sec_mass)
    print('Primary starmass:', pri_starmass, 'Secondary starmass:', sec_SM)
    print('Delta M:', 'g:', delta_M, 'delta DM:', delta_DM, 'delta star:', delta_star, 'delta mass:', delta_mass )
   
    
    return None

def LGal_GLF(index, bin_num):
    #Evaluate the galaxy group paramters with r2t, halo condition limit

    #Import data
    basePath = './TNG100' #set basepath
    path = './TNG100/L_Galaxies/LGalaxies_099.hdf5'
    #Set fields for data import
    FOF_fields = ['GroupFirstSub', 'GroupNsubs', 'GroupPos'] #For FOF halos
    sub_fields = [ 'SubhaloGrNr', 'SubhaloFlag' ,'SubhaloCM', 'SubhaloMass','SubhaloMassType', 'SubhaloVelDisp', 'SubhaloStellarPhotometrics'] #For subfind subhalos
    #data import from TNG100 snapshot 99 (z=0) groupcat
    halos = il.groupcat.loadHalos(basePath, 99, fields= FOF_fields)
    subhalos = il.groupcat.loadSubhalos(basePath,99,fields=sub_fields)    
    L_Gal=  h5py.File(path, 'r')['Galaxy']

    MagDust = np.array(L_Gal['MagDust'])
    TNG_index = np.array(L_Gal['SubhaloIndex_TNG'])

    CentDist = np.array(L_Gal['DistanceToCentralGal'])
    StellarMass = np.array(L_Gal['StellarMass'])
    vel_dis = subhalos['SubhaloVelDisp']
    R2t = vel_dis/368*1000
    #Test for condtions
    
    
    pri = TNG_index[index]
    FOF  = subhalos['SubhaloGrNr'][pri]
    sub_count = halos['GroupNsubs'][FOF]
    pri_mag = MagDust[index]
    pri_mass = subhalos['SubhaloMass'][pri]
    pri_DMmass = subhalos['SubhaloMassType'][pri][1]


    
    pri_R2t = R2t[pri]
    Mag_list = [pri_mag[1]]
    DM_list =[pri_DMmass]

    count = 0

    for sub in range(pri+1, pri+ sub_count):

        if sub in TNG_index:
            L_ind  = np.where(TNG_index == sub)[0][0]

            sub_mag = MagDust[L_ind]
            
            dist  = np.sqrt(CentDist[L_ind][0]**2+CentDist[L_ind][1]**2+CentDist[L_ind][2]**2)
            if dist<= pri_R2t and sub_mag[16]<99:
                count += 1
                sub_mass = subhalos['SubhaloMass'][sub]
                sub_DMmass = subhalos['SubhaloMassType'][sub][1]

                
                DM_list.append(sub_DMmass)
                Mag_list.append(sub_mag[1])
        
  
    
    plt.figure(0)
    fig, axs = plt.subplots(2)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
    fig.suptitle('GLF and GMF')
    axs[0].hist(Mag_list, bins = bin_num)

    axs[0].set_xlabel('$ M_r $')
    axs[0].set_ylabel('Frequency')
    
    axs[1].hist(np.log10(DM_list), bins = bin_num)

    axs[1].set_xlabel('$log DM Mass(log_{10}(M_{sol}*10^{10}))$')
    axs[1].set_ylabel('Frequency')
    plt.gca().invert_xaxis()
    
 
    return None
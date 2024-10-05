# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import requests
import h5py
import matplotlib.pyplot as plt
import numpy as np
baseUrl = 'http://www.tng-project.org/api/'
headers = {'api-key':'e230a557d3f126462f10a705275b9f53'}

def get(path, params=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically
    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string
    return r

API = get(baseUrl)
names = [sim['name'] for sim in API['simulations']]
sim = get( API['simulations'][4]['url'] )
snaps = get( sim['snapshots'] )
snap = get( snaps[-1]['url'] )
subs = get( snap['subhalos'], {'limit':20, 'order_by':'-mass_stars'} )
sub = get( subs['results'][1]['url'] )
url = sub['related']['parent_halo'] + "info.json"
parent_fof = get(url)
mpb1 = get( sub['trees']['sublink_mpb'] ) # file saved, mpb1 contains the filename
f = h5py.File(mpb1,'r')
mpb2 = get( sub['trees']['lhalotree_mpb'] )
with h5py.File(mpb2,'r') as f:
     pos = f['SubhaloPos'][:]
     snapnum = f['SnapNum'][:]
     subid = f['SubhaloNumber'][:]
    
for i in range(3):
     plt.plot(snapnum,pos[:,i] - pos[0,i], label=['x','y','z'][i])
plt.legend()
plt.xlabel('Snapshot Number')
plt.ylabel('Pos$_{x,y,z}$(z) - Pos(z=0)');
plt.close()
url2 = sim['snapshots'] + "z=1/"
snap = get(url2)
i = np.where(snapnum == 85)
sub_prog_url = "http://www.tng-project.org/api/Illustris-3/snapshots/85/subhalos/185/"
sub_prog = get(sub_prog_url)

cutout = get(sub_prog_url+"cutout.hdf5", {'gas':'Coordinates,Masses'})

with h5py.File(cutout,'r') as f:
     x = f['PartType0']['Coordinates'][:,0] - sub_prog['pos_x']
     y = f['PartType0']['Coordinates'][:,1] - sub_prog['pos_y']
     dens = np.log10(f['PartType0']['Masses'][:])
 
plt.hist2d(x,y,weights=dens,bins=[150,100])
plt.xlabel('$\Delta x$ [ckpc/h]')
plt.ylabel('$\Delta y$ [ckpc/h]');
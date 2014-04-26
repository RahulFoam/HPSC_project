# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 21:34:52 2014

@author: numguy
"""
import numpy as np

def weightx(x_w,scheme):
    (r,c) = np.shape(x_w)
    wpx1 = np.zeros((r,c-3))
    wpx2 = np.zeros((r,c-3))
    wpx3 = np.zeros((r,c-3))
    wnx1 = np.zeros((r,c-3))
    wnx2 = np.zeros((r,c-3))
    wnx3 = np.zeros((r,c-3))
    if scheme == 1:
        wpx2 = wpx2 + 1
        wnx2 = wnx2 + 1
    elif scheme == 2:
        wpx2 = (2*x_w[:,1:-2] + x_w[:,0:-3])/(x_w[:,1:-2] + x_w[:,0:-3])
        wnx2 = (2*x_w[:,2:-1] + x_w[:,3:])/(x_w[:,2:-1] + x_w[:,3:])
        wpx3 = (-1*x_w[:,1:-2])/(x_w[:,1:-2] + x_w[:,0:-3])
        wnx3 = (-1*x_w[:,2:-1])/(x_w[:,2:-1] + x_w[:,3:])
    elif scheme == 3:
        wpx1 = (x_w[:,1:-2]*(2*x_w[:,1:-2] + x_w[:,0:-3]))/((x_w[:,1:-2] + x_w[:,2:-1])*(x_w[:,2:-1] + 2*x_w[:,1:-2] + x_w[:,0:-3]))
        wnx1 = (x_w[:,2:-1]*(2*x_w[:,2:-1] + x_w[:,3:]))/((x_w[:,1:-2] + x_w[:,2:-1])*(x_w[:,1:-2] + 2*x_w[:,2:-1] + x_w[:,3:]))
        wpx2 = (x_w[:,2:-1]*(2*x_w[:,1:-2] + x_w[:,0:-3]))/((x_w[:,1:-2] + x_w[:,2:-1])*(x_w[:,1:-2] + x_w[:,0:-3]))
        wnx2 = (x_w[:,1:-2]*(2*x_w[:,2:-1] + x_w[:,3:]))/((x_w[:,1:-2] + x_w[:,2:-1])*(x_w[:,1:-2] + x_w[:,3:]))
        wpx3 = (-1*x_w[:,1:-2]*x_w[:,2:-1])/((x_w[:,1:-2] + x_w[:,0:-3])*(2*x_w[:,2:-1] + x_w[:,1:-2] + x_w[:,0:-3]))
        wnx3 = (-1*x_w[:,1:-2]*x_w[:,2:-1])/((x_w[:,3:] + x_w[:,2:-1])*(x_w[:,1:-2] + 2*x_w[:,2:-1] + x_w[:,3:]))
    return wpx2
 

def weighty(y_w,scheme):
    (r,c) = np.shape(y_w)
    wpy1 = np.zeros((r-3,c))
    wpy2 = np.zeros((r-3,c))
    wpy3 = np.zeros((r-3,c))
    wny1 = np.zeros((r-3,c))
    wny2 = np.zeros((r-3,c))
    wny3 = np.zeros((r-3,c))
    if scheme == 1:
        wpy2 = wpy2 + 1
        wny2 = wny2 + 1
    elif scheme == 2:
        wpy2 = (2*y_w[2:-1,:] + y_w[3:,:])/(y_w[2:-1,:] + y_w[3:,:])
        wny2 = (2*y_w[1:-2,:] + y_w[0:-3,:])/(y_w[1:-2,:] + y_w[0:-3,:])
        wpy3 = (-1*y_w[2:-1,:])/(y_w[2:-1,:] + y_w[3:,:])
        wny3 = (-1*y_w[1:-2,:])/(y_w[1:-2,:] + y_w[0:-3,:])
    elif scheme == 3:
        wpy1 = (y_w[2:-1,:]*(2*y_w[2:-1,:] + y_w[3:,:]))/((y_w[2:-1,:] + y_w[1:-2,:])*(y_w[1:-2,:] + 2*y_w[2:-1,:] + y_w[3:,:]))
        wny1 = (y_w[1:-2,:]*(2*y_w[1:-2,:] + y_w[0:-3,:]))/((y_w[2:-1,:] + y_w[1:-2,:])*(2*y_w[1:-2,:] + y_w[2:-1,:] + y_w[0:-3,:]))
        wpy2 = (y_w[1:-2,:]*(2*y_w[2:-1,:] + y_w[3:,:]))/((y_w[2:-1,:] + y_w[1:-2,:])*(y_w[2:-1,:] + y_w[3:,:]))
        wny2 = (y_w[2:-1,:]*(2*y_w[1:-2,:] + y_w[0:-3,:]))/((y_w[2:-1,:] + y_w[1:-2,:])*(y_w[2:-1,:] + y_w[0:-3,:]))
        wpy3 = (-1*y_w[2:-1,:]*y_w[1:-2,:])/((y_w[2:-1,:] + y_w[3:,:])*(2*y_w[1:-2,:] + y_w[2:-1,:] + y_w[3:,:]))
        wny3 = (-1*y_w[2:-1,:]*y_w[1:-2,:])/((y_w[0:-3,:] + y_w[1:-2,:])*(2*y_w[1:-2,:] + y_w[2:-1,:] + y_w[0:-3,:]))
    return wpy2

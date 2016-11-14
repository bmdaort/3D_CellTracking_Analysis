import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def get_dR_Polarization(xys_s,dxys):
    xyr=[]
    for i in list(range(len(xys_s))):
        xyr.append(xys_s[0]-np.mean(xys_s[0]))
    xyr_r=[]
    dxyr=[]
    thr=[]
    dR=[]
    ang=[]
    for i in list(range(len(dxys))):
        u,s,rm=np.linalg.svd(dxys[i], full_matrices=True)
        rm_a=np.array(rm)
        rm_a=np.transpose(rm_a)
        dxyr.append(np.dot(np.array(dxys[i]),rm_a))
        xyr_r.append(np.dot(np.array(xyr[i]),rm_a))
        thr.append(np.divide(np.arctan2(dxyr[i][:,1],dxyr[i][:,0]),np.pi/180))
        thr[i][thr[i]<0]=thr[i][thr[i]<0]+360
        dR.append(np.sqrt(np.square(dxyr[i][:,0])+np.square(dxyr[i][:,1])))
        ang.append(np.arctan2(dxyr[i][:,1],dxyr[i][:,0]))
    
    dR2=[item for sublist in dR for item in sublist]
    ang2=[item for sublist in ang for item in sublist]
    ang3=pd.Series(ang2)
    dR3=pd.Series(dR2)
    angpoint=np.linspace(-np.pi,np.pi,25);
    angstep=angpoint[1::]/2+angpoint[0:-1]/2;
    rho=[]
    dR_thetha=[]
    for i in list(range(len(angpoint)-1)):
        rho.append(np.mean(dR3.loc[np.logical_and(ang3>=angpoint[i],ang3<angpoint[i+1])]))
        dR_thetha.append(angstep[i])
    axx=rho*np.cos(dR_thetha)/np.mean(rho)
    ayy=rho*np.sin(dR_thetha)/np.mean(rho)
    #return dR_thetha
    return {'axx':axx, 'ayy':ayy}
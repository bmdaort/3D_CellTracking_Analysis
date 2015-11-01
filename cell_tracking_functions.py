import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def get_xys(df):
    N=df.Nc
    Nc=N.loc[len(N)-1]#Number of cells
    xys=[]
    traj_len=[] #lenght of each trajectory
    for i in list(range(1,Nc+1)):
        c=df[df['Nc']==i]
        xys.append(c[['X','Y']])
        traj_len.append(len(c))
    Nmin=min(traj_len)
    return {'xys':xys, 'Nmin':Nmin, 'Nc':Nc}

def get_xys_s(xys,Nmin):
    xys_s=[] #xys_s is just xys sub-sampled
    for i in list(range(len(xys))):
        xys_s.append(xys[i][0:Nmin])
    return xys_s
    
def get_timeVector(df,Nmin):
    tframe=df['tframe']
    t_min=[]
    dt=tframe[0]/3600#Time in hours
    for i in list(range(Nmin-1)):
        t_min.append(dt*i)
    t_min=pd.Series(t_min)
    return t_min

def get_dxys(xys_s):
    dxys=[]
    for i in list(range(len(xys_s))):
        dxys.append(np.array(xys_s[i].iloc[1:,0:])-np.array(xys_s[i].iloc[0:-1,0:]))
    return dxys

def get_velocity(dxys):
    v=[]
    v2=[]
    for i in list(range(len(dxys))):
    	vels=np.sqrt(np.square(dxys[i][:,0])+np.square(dxys[i][:,1]))
        v.append(vels)
        v2.append(np.mean(vels))
    v_mean=np.mean(v,axis=0,dtype=np.float32)
    v_std=np.std(v,axis=0,dtype=np.float32)
    v_mean=pd.Series(v_mean)
    v_std=pd.Series(v_std)
    return {'v_Mean':v_mean, 'v_Std':v_std, 'velocities':v2}

def sub_sample_vel(v_mean,t_min,window):
    v_mean_window=[]
    t_min_window=[]
    for i in list(range(len(v_mean)/window)):
        v_mean_window.append(np.mean(v_mean[i*(window):i*(window)+(window-1)]))
        t_min_window.append(t_min[i*(window)+(window-1)])
    return {'v_Mean_subsample':v_mean_window, 't_min_subsample':t_min_window}

def get_acf(dxys):
    M=len(dxys[0])
    N=len(dxys)
    acf=[]
    for j in list(range(N)):
        c_acf=[]
        for i in list(range(M-1)):
            c_acf.append(np.divide(sum(sum(dxys[j][0:M-i]*dxys[j][i:M])),(M-i)))
        acf.append(c_acf)
    acf_mean=np.mean(acf,axis=0,dtype=np.float32)
    acf_mean=acf_mean/acf_mean[2];
    acf_mean=acf_mean[2:100]
    return acf_mean

def get_msd(xys_s):
    M=len(xys_s[0])
    N=len(xys_s)
    msd=[]
    for j in list(range(N)):
        msdr=[]
        for i in list(range(M-1)):
            msd1=np.array(xys_s[j].iloc[i:M,:])-np.array(xys_s[j].iloc[0:M-i,:])
            msd2=np.mean(np.square(msd1),axis=0,dtype=np.float32)
            msd3=sum(msd2)
            msdr.append(msd3)
        msd.append(msdr)
    msd_mean=np.mean(msd,axis=0,dtype=np.float32)
    return {'Mean':msd_mean, 'Raw':msd}

def get_dR_PDF(dxys):
    dxmax=50
    binn=70
    bins=np.linspace((dxmax/binn)/2, dxmax, binn)
    feq, binw = np.histogram(np.absolute(dxys), bins,density=True)
    return {'feqw':feq, 'bin':binw}
	
def get_dThetha_PDF(xys_s):
    tloi=[1,5,10]
    dth000=[]
    for i in list(range(len(xys_s))):
        dth00=[]
        for j in list(range(len(tloi))):
            c_xys=[]
            dxyt=[]
            dr=[]
            dth=[]
            dth0=[]
            c_xys.append(xys_s[i][0::tloi[j]])
            dxyt.append(np.array(c_xys[0].iloc[1:,0:])-np.array(c_xys[0].iloc[0:-1,0:]))
            dr.append(np.sqrt(np.sum(np.square(dxyt[0]),axis=1)))
            dth.append(np.divide(np.sum(np.multiply(dxyt[0][1::],dxyt[0][0:-1:]),axis=1),np.multiply(dr[0][1::],dr[0][0:-1:])))
            dth0.append(np.arccos(dth[0]))
            dth00.append(np.divide(dth0[0],np.pi/180))
        dth000.append(dth00)
    dthetha=[]
    for j in list(range(len(tloi))):
        dthetha_l=[]
        for i in list(range(len(xys_s))):
            dthetha_l.append(dth000[i][j].tolist())
        dthetha.append([item for sublist in dthetha_l for item in sublist])
        
    n_deg_bin=np.linspace((180/6)/2, 180, 6) # number of angles to consider for histogram
    dthetha_dist=[]
    for i in list(range(len(tloi))):
        dthetha_array=np.array(dthetha[i])
        dthetha_array=dthetha_array[~np.isnan(dthetha_array)]
        feq_thetha, bin_thetha=np.histogram(dthetha_array,5,density=True)
        dthetha_dist.append(feq_thetha)    
        
    return {'thetha_dist':dthetha_dist, 'bin_theta':bin_thetha}







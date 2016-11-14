import cell_tracking_functions as ct
import numpy as np
import matplotlib as plt

def get_Velocities_for_DataFrame(DF):
    GetXYS=ct.get_xys(DF)
    GetVel=ct.get_velocity(ct.get_dxys(GetXYS['xys']))#Get velocity
    return GetVel

def get_Velocities_for_DataFrame_noMeans(DF):
    GetXYS=ct.get_xys(DF)
    GetVelnoMeans=ct.get_velocity_noMeans(ct.get_dxys(GetXYS['xys']))#Get velocity
    return GetVelnoMeans

def get_Velocities_for_Dataframe_at_lag(DF):
    dt=DF['tframe'][0]/60
    tlag_n=int(tlag/dt)
    GetXYS=ct.get_xys(DataFrames[i][j]) #This returns the XY trajectories for the specific DataFrame
    xys_s_temp=ct.get_xys_s(GetXYS['xys'],GetXYS['Nmin'])
    GetVel_lag=ct.get_velocity_lag(xys_s_temp,tlag_n)
    return GetVel_lag

def get_Total_Invdist_for_DataFrame(DF):
    GetXYS=ct.get_xys(DF)
    invdist=ct.get_invdist(GetXYS['xys'])
    return invdist

def get_Path_length_for_DataFrame(DF):
    GetXYS=ct.get_xys(DF)
    CumulativeDistance=ct.get_invdist_Cummulative(ct.get_xys_s(GetXYS['xys'],GetXYS['Nmin']))
    return CumulativeDistance

def get_Max_Inv_dist_for_DataFrame(DF):
    GetXYS=ct.get_xys(DF)
    MaxInvDist=ct.get_Max_Inv_dist(GetXYS['xys'])
    return MaxInvDist

def get_trajectories_singleAxis_for_DF(DF):
    from random import randint
    sns.set_palette(sns.color_palette("Paired"))
    fig = plt.figure(figsize=(5, 5),frameon=False)
    ax=fig.add_subplot(1, 1, 1)
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    #ax.spines['left'].set_smart_bounds(True)
    #ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylim([-300,300])
    ax.set_xlim([-300,300])
    ticklab = ax.xaxis.get_ticklabels()[0]
    ax.xaxis.set_label_coords(300, -40,transform=ticklab.get_transform())
    ax.set_xlabel('x($\mu$m)',fontsize=14)
    ticklab = ax.yaxis.get_ticklabels()[0]
    ax.yaxis.set_label_coords(90, 280,transform=ticklab.get_transform())
    ax.set_ylabel('y($\mu$m)',rotation=0,fontsize=14)
    
    
    GetXYS=ct.get_xys(DF)
    xys_s=ct.get_xys_s(GetXYS['xys'],GetXYS['Nmin'])
    for n in list(range(12)):
        m=randint(0,len(xys_s))-1
        xys_s_x_n=xys_s[m]['X']-(xys_s[m]['X'][xys_s[m]['X'].index[0]])
        xys_s_y_n=xys_s[m]['Y']-(xys_s[m]['Y'][xys_s[m]['X'].index[0]])
        xys_s_x_n=[x*(100/(60.)) for x in xys_s_x_n]
        xys_s_y_n=[x*(100/(60.)) for x in xys_s_y_n]
        ax.plot(xys_s_x_n,xys_s_y_n)
    return fig

def get_trajectories_for_DF(DF):
    GetXYS=ct.get_xys(DF)
    xys_s=ct.get_xys_s(GetXYS['xys'],GetXYS['Nmin'])
    plt.figure(figsize=(5, 5),frameon=False)
    for m in list(range(9)):
        plt.plot()
        plt.subplot(3,3,m+1)
        xys_s_x_n=xys_s[m]['X']-min(xys_s[m]['X'])
        xys_s_y_n=xys_s[m]['Y']-min(xys_s[m]['Y'])
        plt.plot(xys_s_x_n,xys_s_y_n)
        plt.axis('off')
        axes = plt.gca()
        axes.set_ylim([0,125])
        axes.set_xlim([0,125])

def get_MSD_for_DataFrame(DF):
    GetXYS=ct.get_xys(DF)
    xys_s_temp=ct.get_xys_s(GetXYS['xys'],GetXYS['Nmin'])
    i=0
    for j in list(range(len(xys_s_temp))):
        if np.isnan(xys_s_temp[i]['X']).any():
            xys_s_temp.pop(i)
            i=i-1
        i=i+1
    msd=ct.get_msd(xys_s_temp)['Mean']#See get_msd function
    return msd

def get_MSD_for_DataFrame_raw(DF):
    GetXYS=ct.get_xys(DF)
    xys_s_temp=ct.get_xys_s(GetXYS['xys'],GetXYS['Nmin'])
    i=0
    for j in list(range(len(xys_s_temp))):
        if np.isnan(xys_s_temp[i]['X']).any():
            xys_s_temp.pop(i)
            i=i-1
        i=i+1
    msd=ct.get_msd(xys_s_temp)['Raw']#See get_msd function
    return msd

def get_MSD_raw(DF):
    GetXYS=ct.get_xys(DF)
    xys_s_temp=GetXYS['xys']
    i=0
    for j in list(range(len(xys_s_temp))):
        if np.isnan(xys_s_temp[i]['X']).any():
            xys_s_temp.pop(i)
            i=i-1
        i=i+1
    msd=ct.get_msd(xys_s_temp)['Raw']#See get_msd function
    return msd


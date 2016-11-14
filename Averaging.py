import numpy as np
def get_AVG(n):
    ld=[]
    for k in list(range(len(n))):
        ld.append(len(n[k]))
    l_min=min(ld)
    for k in list(range(len(n))):
        n[k]=n[k][0:l_min]
    n_avg=np.mean(n,axis=0)
    n_std=np.std(n,axis=0)
    return {'Mean':n_avg, 'STD':n_std, 'l_min':l_min}
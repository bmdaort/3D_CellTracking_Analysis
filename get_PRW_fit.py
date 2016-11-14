import numpy as np
import pandas as pd
import scipy
import lmfit as lm

def PRW(x, S,P,d):
	return np.square(S) *2*P* (x - P * ( 1 - np.exp(-x/P) ) ) + 4 * d
	
def get_PRW_fit(msd,t_min):
	t=list(range(1,len(t_min)/3))
	t=np.array(t,dtype=float)
	Nt=len(t);
	wif=((2*np.square(t)+1)/3)/t/(Nt-t+1)
	r_squared=[]
	for i in list(range(len(msd))):
		wt=np.sqrt(1/np.square(wif)/msd[i][1:len(t_min)/3])
		time = t_min[1:len(t_min)/3]
		y = msd[i][1:len(t_min)/3]
		gmod = lm.Model(PRW)
		gmod.set_param_hint('S', value = 5, min=0, max=1000)
		gmod.set_param_hint('P', value = 1, min=0, max=20)
		gmod.set_param_hint('d', value = 1, min=0, max=100)
		result = gmod.fit(y,x=time,weights=wt)
		r_squared.append(1-result.residual.var()/np.var(y))
	return r_squared
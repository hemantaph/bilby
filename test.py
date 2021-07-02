import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 

path="hplog.csv"
df=pd.read_csv(path)
f= np.asarray(df['Frequency'].values)
hp= np.asarray(df['hp'].values)

f = f.astype(np.float)
hp = hp.astype(np.float)

plt.loglog(f,hp)

plt.savefig('img.png')
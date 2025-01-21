import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('./outputs/output.txt', names=['ID','x','y'])

p1 = data[data['ID'] == 0]
p2 = data[data['ID'] == 1]

slice_int = 500

plt.plot(p1.index/2, p2['y'].to_numpy()-p1['y'].to_numpy())
plt.plot(p1.index/2, 16 * np.ones_like(p1['y']))
# plt.plot(p1.index, p1['y'])
# plt.plot(p2.index, p2['y'])
# plt.plot(p1.index[:slice_int], p1['y'][:slice_int])
# plt.plot(p2.index[:slice_int], p2['y'][:slice_int])
plt.show()

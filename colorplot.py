import numpy as np
import matplotlib.pyplot as plt

P = [0,1.1,1.4]
chi=np.zeros([100,10])
for i in range(100):
    for j in range(100):
        z[i,j] = i+j

plt.imshow(chi)
plt.colorbar () # カラーバーの表示 
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

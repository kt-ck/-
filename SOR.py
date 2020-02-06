# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 19:10:02 2019
主题：SOR
内容: w取0.3、0.6、0.9、1.2、1.5、1.8
@author: KT
"""

import numpy as np
#方程组

A = np.array([[4, -1, 0, -1, 0, 0],
              [-1, 4, -1, 0,-1, 0],
              [0, -1,  4, 0, 0,-1],
			  [-1, 0,  0, 4, -1,0],
              [0, -1,  0,-1, 4,-1],
			  [0, 0,  -1, 0,-1, 4]])

b = np.array([2,3,2,2,1,2])

L = np.zeros((6, 6))

U = np.zeros((6, 6))

D = np.eye(6)

X_pre = np.zeros(6)

X_cur = np.zeros(6)

#计算 L
for i in range(1, 6):
    for j in range(0,i):
        L[i, j] = A[i, j] * -1
        

#计算 U
for i in range(0, 5):
    for j in range(i + 1, 6):
        U[i, j] = A[i, j] * -1
      
#计算 D
for i in range(0, 6):
    D[i, i] = A[i, i]; 

w = 1.110
x = []
y = []
while 1.110 <= w <= 1.2:
    #计算 SOR迭代矩阵
    J = np.linalg.inv(D-w * L).dot((1 - w) * D + w * U)
        
    #计算 SOR迭代向量
    f = np.linalg.inv(D - w * L).dot(b) * w
        
    #SOR迭代
    #print("\n\nw = " + str(w) + '\n')
    X_cur = J.dot(X_pre) + f
        
    cnt = 1
     
    while np.linalg.norm(X_cur - X_pre, ord = np.inf) > 10E-5:
        X_pre = X_cur
        X_cur = J.dot(X_pre) + f
        cnt = cnt + 1
        
    #print("SOR迭代一共做了"+str(cnt)+"次，最终的近似解为：")       
    #print(X_cur)
    x.append(w)
    y.append(cnt)    
    w = w + 0.001
    X_pre = np.zeros(6)
    X_cur = np.zeros(6)    
print(x)
print(y)

import matplotlib.pyplot as plt

plt.plot(x,y)
plt.show()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

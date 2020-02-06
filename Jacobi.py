# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 13:25:52 2019
主题：Jacobi迭代
@author: CK
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

#求迭代矩阵
J = np.linalg.inv(D).dot(L + U)

#jacobi迭代向量
f = np.linalg.inv(D).dot(b)

print("Jacobi迭代矩阵为：")
print(J)
print("\nJacobi迭代向量为：")
print(f)
#进行jacobi迭代
print("\n迭代过程：")
X_cur = J.dot(X_pre) + f

cnt = 1                #迭代次数

print("\n第"+str(cnt)+"次迭代后的近似解为:")
print(X_cur)
print("此时前后两次迭代矩阵的无穷范数为：")
print(np.linalg.norm(X_cur - X_pre, ord = np.inf))

while np.linalg.norm(X_cur - X_pre, ord = np.inf) > 10E-5:
    X_pre = X_cur
    X_cur = J.dot(X_pre) + f
    cnt = cnt + 1
    print("\n第"+str(cnt)+"次迭代后的近似解为:")
    print(X_cur)
    print("此时前后两次迭代矩阵的无穷范数为：")
    print(np.linalg.norm(X_cur - X_pre, ord = np.inf))
    
print("\n\n迭代结束。")
print("jacobi迭代一共做了"+str(cnt)+"次，最终的近似解为：")
print(X_cur)

























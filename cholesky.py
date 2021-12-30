# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 01:36:51 2021

@author: paes mario
"""

#Let be the matrix A, calculate the Cholesky decomposition by hand
A = [[ 16, 4, 8, 4],
     [ 4, 10, 8, 4],
     [ 8, 8, 12, 10],
     [ 4, 4, 10, 12]]


'''First we need to verfy if A matrix is positive definite
and I will use the Sylvester's criterion
Sylvester's criterion states that a symmetric
(more generally, Hermitian) matrix is positive definite
if and only if its principal minors are all positive
So I will check by hand if A matches the criterion
If not, the program will return an alert'''

#writing the functions for Sylvester's criterion
def det_1(x):
    return x[0][0]

def det_2(x):
    a=x[0][0]*x[1][1]
    b=x[0][1]*x[1][0]
    return a-b

def det_3(x):  
    x.append(x[0]);x.append(x[1]);a=0    
    for i in range(3):        
        r=1
        for j in range(3):            
            r*=x[i+j][j]
        a+=r
    b = 0
    for k in range(3):        
        s = 1
        aux = 0
        for l in range(2,-1,-1):            
            s*= x[k+aux][l]            
            aux+=1
        aux+=1;b+=s
    return a-b

def Ak(x):
    xi,xm,xf=[],[],[]
    for i in range(len(x)-1):
        xi.append(x[i])
    for j in range(len(x)-1):
        xm=[]
        for k in range(len(x)-1):
            xm.append(xi[j][k])
        xf.append(xm)
    return xf

def all_Ak_pos(x):
    x3=Ak(x)
    x2=Ak(x3)
    x1=Ak(x2)
    
    return (det_1(x1) and det_2(x2) and det_3(x3))>0

def simetric(x):
    aux=0
    for i in range(len(x)):
        for j in range(len(x)):
            if x[i][j]!=x[j][i]:
                aux+=1
    return aux==0

def Sylvester(x):
    
    return all_Ak_pos(x) and simetric(x)

'''Let's create a Cholesky decomposition by hand
for this purpose
we need a couple of helping functions'''

#absolute number
def abs(x):
    if x>0:
        abs=x
    else:
        abs=-x
    return abs

#squared root
def sqrt(x):
    root = x
    tol = 1e-5
    while abs(x - root * root) > tol:
        root = (root + x / root) / 2
    return root

#zero matrix nxn
def mx_zero(x):
    mx=[]
    for i in range(x):
        aux=[]
        for j in range(x):
            aux.append(0)
        mx.append(aux)
    return mx

'''Now let's use this helping functions
to code my Cholesky decomposition function'''

#Cholesky decomposition
def Cholesky(x):
    G=mx_zero(len(x))
    for i in range(len(x)):
        for j in range(i+1):
            add = sum(G[i][k] * G[j][k] for k in range(j))
            if (i==j):
                G[i][j]=sqrt(x[i][i]-add)
            else:
                G[i][j]=(1.0 / G[j][j] * (x[i][j]-add))
    return G

'''Now we calculate the Cholesky decomposition of matrix A
And include Sylvester criterion'''

if Sylvester(A):
    Result = Cholesky(A)
    print("My Cholesky decomposition by hand is: \n", Result)
    #Validating the result with Numpy Library
    
    import numpy as np
    
    Arr = np.array(A)
    Answer = np.linalg.cholesky(A)
    
    print("Numpy answer is: \n", Answer)
    
else:
    
    print("Unable to calculate the Cholesky decomposition \n"+
"Please, check if the matrix A is positive definite!")

'''
Output:
    
My Cholesky decomposition by hand is: 
 [[4.000000636692939, 0, 0, 0], 
  [0.9999998408267905, 3.0000000544547163, 0, 0], 
  [1.999999681653581, 2.0000001759277817, 2.0000002353410418, 0], 
  [0.9999998408267905, 1.0000000879638908, 2.9999997894070662, 1.000001405976164]]
Numpy answer is: 
 [[4. 0. 0. 0.]
 [1. 3. 0. 0.]
 [2. 2. 2. 0.]
 [1. 1. 3. 1.]]
'''

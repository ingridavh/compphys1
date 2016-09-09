#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

myfile1 = open('../build-FYS4150_project1-Desktop-Debug/p1_result_N10.txt')
myfile2 = open('../build-FYS4150_project1-Desktop-Debug/p1_result_N100.txt')
myfile3 = open('../build-FYS4150_project1-Desktop-Debug/p1_result_N1000.txt')

N1 = 10
N2 = 100
N3 = 1000

numbers1 = []
numbers2 = []
numbers3 = []

x1 = np.linspace(0, 1, N1+2)
x2 = np.linspace(0, 1, N2+2)
x3 = np.linspace(0, 1, N3+2)

v1 = np.zeros(N1+2)
v2 = np.zeros(N2+2)
v3 = np.zeros(N3+2)

v_ex1 = np.zeros(N1+2)
v_ex2 = np.zeros(N2+2)
v_ex3 = np.zeros(N3+2)

for line in myfile1:
    line = line.split()
    if line:
        line = [float(i) for i in line]
        numbers1.append(line)

for line in myfile2:
    line = line.split()
    if line:
        line = [float(i) for i in line]
        numbers2.append(line)

for line in myfile3:
    line = line.split()
    if line:
        line = [float(i) for i in line]
        numbers3.append(line)


for j in range(N1+2):
    row = numbers1[j]
    v1[j] = row[0]
    v_ex1[j] = row[1]

for j in range(N2+2):
    row = numbers2[j]
    v2[j] = row[0]
    v_ex2[j] = row[1]

for j in range(N3+2):
    row = numbers3[j]
    v3[j] = row[0]
    v_ex3[j] = row[1]


plt.subplot(3, 1, 1)
plt.plot(x1, v1, '-r', x1, v_ex1, '-b')
plt.title('N=10')
plt.xlabel('x')
plt.ylabel('v(x)')

plt.subplot(3, 1, 2)
plt.plot(x2, v2, '-r', x2, v_ex2, '-b')
plt.title('N=100')
plt.xlabel('x')
plt.ylabel('v(x)')

plt.subplot(3, 1, 3)
plt.plot(x3, v3, '-r', x3, v_ex3, '-b')
plt.title('N=1000')
plt.xlabel('x')
plt.ylabel('v(x)')

plt.show()

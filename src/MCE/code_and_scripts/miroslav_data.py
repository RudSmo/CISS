#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 12:14:35 2020

@author: hopjan
"""
import numpy as np
import scipy.linalg as sl
import scipy.special as sp
import matplotlib.pyplot as plt

########################## data ED code ########################################
###############################################################################
###############################################################################
###############################################################################
time=np.arange(0, 120, 0.1)
#.*L2M40t01.0omega0.*gam.*dt0.10tot.*PB0.bin

gamma = 0.4		# 0.05 0.1, 0.2, 0.4
n0=np.fromfile("dimer_data/n0L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)
n1=np.fromfile("dimer_data/n1L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)
x0=np.fromfile("dimer_data/x0L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)
p0=np.fromfile("dimer_data/p0L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)
x1=np.fromfile("dimer_data/x1L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)
p1=np.fromfile("dimer_data/p1L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)

######################### density plots #######################################
###############################################################################
###############################################################################
###############################################################################

plt.figure(1)
plt.plot(time, n1, linestyle='-')
plt.xlabel('t')
plt.ylabel('n_1')


plt.figure(2)
plt.plot(time, n0, linestyle='-')
plt.xlabel('t')
plt.ylabel('n_2')


plt.figure(10)
plt.plot(time, n0, label="site 0")
plt.plot(time, n1, label="site 1")
plt.xlabel('t')
plt.ylabel('n')
plt.legend()


########################## p, p_q plots #######################################
###############################################################################
###############################################################################
###############################################################################

plt.figure(3)
plt.plot(time, x1-x0, linestyle='-')
plt.xlabel('t')
plt.ylabel('q')

plt.figure(4)
plt.plot(time, p1-p0, linestyle='-')
plt.xlabel('t')
plt.ylabel('p_q')

plt.figure(5)
plt.plot(-x1+x0, -p1+p0, linestyle='-')
plt.xlabel('q')
plt.ylabel('p_q')

plt.figure(6)
plt.plot(-x0, -p0, linestyle='-')
plt.xlabel('x_1')
plt.ylabel('p_1')

plt.figure(7)
plt.plot(-x1, -p1, linestyle='-')
plt.xlabel('x_2')
plt.ylabel('p_2')

########################## data ED code ########################################
###############################################################################
###############################################################################
###############################################################################

n0=np.fromfile("dimer_data_polstate/n0L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)
n1=np.fromfile("dimer_data_polstate/n1L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)
x0=np.fromfile("dimer_data_polstate/x0L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)
p0=np.fromfile("dimer_data_polstate/p0L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)
x1=np.fromfile("dimer_data_polstate/x1L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)
p1=np.fromfile("dimer_data_polstate/p1L2M40t01.0omega0.1gam%.2fdt0.10tot120PB0.bin"%gamma, dtype=float)


######################### density plots #######################################
###############################################################################
###############################################################################
###############################################################################

plt.figure(1)
plt.plot(time, n1, linestyle='-')
plt.xlabel('t')
plt.ylabel('n_1')


plt.figure(2)
plt.plot(time, n0, linestyle='-')
plt.xlabel('t')
plt.ylabel('n_2')


########################## p, p_q plots #######################################
###############################################################################
###############################################################################
###############################################################################

plt.figure(3)
plt.plot(time, x1-x0, linestyle='-')
plt.xlabel('t')
plt.ylabel('q')

plt.figure(4)
plt.plot(time, p1-p0, linestyle='-')
plt.xlabel('t')
plt.ylabel('p_q')

plt.figure(5)
plt.plot(-x1+x0, -p1+p0, linestyle='-')
plt.xlabel('q')
plt.ylabel('p_q')

plt.figure(6)
plt.plot(-x0, -p0, linestyle='-')
plt.xlabel('x_1')
plt.ylabel('p_1')

plt.figure(7)
plt.plot(-x1, -p1, linestyle='-')
plt.xlabel('x_2')
plt.ylabel('p_2')

plt.show()

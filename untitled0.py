# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:27:34 2022

@author: User
"""
import numpy as np

from projekt1 import *

grs80 = transformacje(model = 'grs80')

wsp = 'wsp_inp.txt'

tablica1 = np.genfromtxt(wsp, delimiter = ',', skip_header = 4)

rzad = 12
kol = 3


blh = np.zeros((rzad, kol))
xy92 = np.zeros((rzad, 2))
xy00 = np.zeros((rzad, 2))
neu = np.zeros((rzad, kol))
katy_odl = np.zeros((rzad, 4))
tablica_all = np.zeros((rzad, 14))

for el in range(rzad):
    blh[el] = grs80.hirvonen(tablica1[el, 0], tablica1[el, 1], tablica1[el, 2])
    xy92[el] = grs80.uklad1992(blh[el, 0], blh[el,1])
    xy00[el] = grs80.uklad2000(blh[el, 0], blh[el,1])
    neu[el] = grs80.xyz2neu(tablica1[el, 0], tablica1[el, 1], tablica1[el, 2], tablica1[el, 0] + 46, tablica1[el, 1] + 10, tablica1[el, 2] + 8)
    katy_odl[el] = grs80.katy_odl(tablica1[el, 0], tablica1[el, 1], tablica1[el, 2], tablica1[el, 0] + 46, tablica1[el, 1] + 10, tablica1[el, 2] + 8)
    
    tablica_all[el] = np.hstack([blh[el], xy92[el], xy00[el], neu[el], katy_odl[el]])


np.savetxt("wszystkie_wsp.txt", tablica_all, delimiter = ',', fmt = ['%.7f', '%12.7f', '%9.3f', '%13.3f', '%13.3f', '%13.3f', '%13.3f', '%13.7f', '%13.7f', '%13.7f', '%13.7f', '%13.7f', '%8.3f', '%8.3f'], encoding = 'UTF', header = 'Zestawienie wszystkich współrzędnych - Agata Wyrzykowska, indeks: 312109 \n-- fi -------- lam ------- h ------- x1992 ------- y1992 ------- x2000 ------- y2000 ---------- n ----------- e ----------- u ---------- Az -------- alfa ------ odl2D -- odl3D--')
import math as m
import numpy as np
import geopandas as geo
import pandas as panda
#from mpl_toolkits import mplot3d
import matplotlib.pyplot as pl
import plotly.express as px

#lot z London Heathrow do Warsaw Frederic Chopin
airport_fi = -0.4605
airport_lambda = 51.4708
airport_h = 25.0


def skos_az_odl(n,e,u):
    tan_A = e/n
    s = m.sqrt(n*n+e*e+u*u)
    cos_z = u/(m.sqrt(n*n+e*e+u*u))
    return tan_A, s, cos_z


def geo2xyz(fi, lam, h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a/m.sqrt(1-e2*np.sin(fi)**2)

    x = (N+h)*np.cos(fi)*np.cos(lam)
    y = (N+h)*np.cos(fi)*np.sin(lam)
    z = (N*(1-e2)+h)*np.sin(fi)
    return x, y, z


def geo2neu(F1, L1, H1, F2, L2, H2):
    a = 6378137  # metry
    e2 = 0.0066943800290  # bez jednostek

    #zamiana wsp lotniska i samolotu na xyz
    pkt_1 = geo2xyz(F1, L1, H1, a, e2)
    pkt_2 = geo2xyz(F2, L2, H2, a, e2)

    #zamiana stopni na radiany
    F1 = np.deg2rad(F1)
    L1 = np.deg2rad(L1)

    R = np.array(
                [[-m.sin(F1) * m.cos(L1), -m.sin(L1), m.cos(F1) * m.cos(L1)],
                 [-m.sin(F1) * m.sin(L1), m.cos(L1), m.cos(F1) * m.sin(L1)],
                 [m.cos(F1), 0, m.sin(F1)]
                 ])

    R =  R.transpose()
    x = np.array([[pkt_2[0] - pkt_1[0]],
                  [pkt_2[1] - pkt_1[1]],
                  [pkt_2[2] - pkt_1[2]]])
    #Mnozenie stransopowanej macierzy R przez x
    neu = R @ x
    return neu


fly = np.loadtxt("lot.txt")

#Zapisanie danych jako array Geopandasowy
data = panda.DataFrame(fly)
air_line = geo.GeoDataFrame(data, geometry=geo.points_from_xy(data[1], data[0]))

#wyswietlenie układu geodezyjnego
fig = px.scatter_geo(air_line, title="Lot z Londynu do Warszawy", lon = air_line.geometry.x,
    lat= air_line.geometry.y)
fig.update_geos(showcountries=True)
fig.show()


#2 część (wykres 3D neu)
three_d = []
n=[]
e=[]
u=[]
#237 = ilosc wierszy air_line
for i in range(0,237):
        three_d.append(geo2neu(airport_fi,airport_lambda,airport_h,air_line[1][i],air_line[0][i],air_line[2][i]))
        n.append(three_d[i][0][0])
        e.append(three_d[i][1][0])
        u.append(three_d[i][2][0])

#wykres 3D neu
fig2 = pl.figure("Wykres neu")
ax = pl.axes(projection='3d')
ax.plot3D(n,e,u,'blue')
ax.scatter3D([n[0],n[48],n[236]],[e[0],e[48],e[236]],[u[0],u[48],u[236]], cmap='Greens')
ax.text(n[0],e[0],10000,"Londyn",None)
ax.text(n[236],e[236],u[236],"Warszawa",None)
ax.text(n[48],e[48],u[48],"Horyzont",None)
ax.set_xlabel('N')
ax.set_ylabel('E')
ax.set_zlabel('U')
ax.ticklabel_format(axis="both", style="sci",scilimits = (0,0))
pl.show()

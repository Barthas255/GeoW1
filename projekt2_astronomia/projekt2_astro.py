from astropy.time import Time
import numpy as np
from astropy.coordinates import Angle
from astropy import units as u
import matplotlib.pyplot as pl
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

#Obserwatorzy
obs1 = (52.232222, 21.008333) # Warszawa fi,lambda
obs2 = (0, 106.828611) # Dżakarta
obs3 = (-35.300000, 149.116666) # Canberra

#Gwiazdozbiór Wodnika, gwiazda Sadalsuud (Beta Aquarii)
rektascencja = 21.525981 # godziny ,alfa
deklinacja = -5.571175 #stopnie, delta

#kat godzinny
def kat_godz(lamb,alfa,h):
    t = Time(['2020-01-28T00:00:00'], format='isot', scale='utc')
    gmst = t.sidereal_time('apparent',longitude='greenwich') # GMST = const, bo dzien sie nie zmienia
    gmst = Angle(gmst, u.rad)
    gmst = float(gmst.to_string(unit=u.degree, decimal=True)) # stopnie
    ut1 = h*1.002737909350795 #godziny
    #obliczenie czasu gwiazdowego(w stopniach)
    S = ut1*15 + lamb + gmst
    #obliczenie kąta godzinnego(w stopniach)
    t = S - alfa * 15
    #print("godzina = ", h,"kat godzinny = ",t, "czas gwiazdowy = ", S)
    return t


#Azymut i zenit
def az_and_zenit(fi, delta, t):
    #zamiana na radiany
    fi = np.deg2rad(fi)
    delta = np.deg2rad(delta)
    t = np.deg2rad(t)
    #zenit
    z = np.arccos(np.sin(fi)*np.sin(delta) + np.cos(fi)*np.cos(delta)*np.cos(t))
    z = np.rad2deg(z)
    #azymut
    licznik = -np.cos(delta)*np.sin(t)
    mianownik = np.cos(fi)*np.sin(delta) - np.sin(fi)*np.cos(delta)*np.cos(t)
    Az = np.arctan(licznik / mianownik)
    Az = np.rad2deg(Az)

    if mianownik < 0: Az += 180
    elif licznik < 0: Az += 360
    if Az > 360:
        Az -= 360
    #print("azymut =",Az,"zenit=", z)
    return z, Az

#Zamiana na xyz
def x_y_z(z, Az, r = 1):
    z = np.deg2rad(z)
    Az = np.deg2rad(Az)
    x = r*np.sin(z)*np.cos(Az)
    y = r * np.sin(z) * np.sin(Az)
    z = r * np.cos(z)
    return x, y, z


'''t = kat_godz(obs1[1],rektascencja,0, 3)
print("kat godzinny: ", t)
xd2 = az_and_zenit(obs1[0],deklinacja,t)
print('odleglosc zenitalna: %f, azymut gwiazdy = %f'%(xd2[0],xd2[1]))
xd3 = x_y_z(xd2[0],xd2[1])
print('X: %f, Y = %f, Z = %f'%(xd3[0], xd3[1], xd3[2]))'''

#przesuwanie utworzonego wykresu
def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)

#Wykresy

x = []
y = []
z = []
x1 = []
y1 = []
z1 = []
x2 = []
y2 = []
z2 = []
h = []
h1 = []
h2 = []

for i in range(1, 25):
    #kat godzinny 3 lokacji
    t = kat_godz(obs1[1],rektascencja, i)
    t1 = kat_godz(obs2[1],rektascencja, i)
    t2 = kat_godz(obs3[1], rektascencja, i)

    #zenit i azymut
    zenitAz = az_and_zenit(obs1[0],deklinacja,t)
    zenitAz1 = az_and_zenit(obs2[0], deklinacja, t)
    zenitAz2 = az_and_zenit(obs3[0], deklinacja, t)

    #wysokość nad horyzontem
    h.append(90 - zenitAz[0])
    h1.append(90 - zenitAz1[0])
    h2.append(90 - zenitAz2[0])

    #wsp. xyz
    xyz = x_y_z(zenitAz[0], zenitAz[1])
    xyz1 = x_y_z(zenitAz1[0], zenitAz1[1])
    xyz2 = x_y_z(zenitAz2[0], zenitAz2[1])

    #xyz w listach
    x.append(xyz[0])
    y.append(xyz[1])
    z.append(xyz[2])
    x1.append(xyz1[0])
    y1.append(xyz1[1])
    z1.append(xyz1[2])
    x2.append(xyz2[0])
    y2.append(xyz2[1])
    z2.append(xyz2[2])

def draw_2D(czas, h, clr, ls_type):
    fig = pl.figure("Wysokość w czasie", figsize=(7,5))
    pl.plot(czas, h, color=clr, ls = ls_type)
    pl.grid(True)
    pl.xlabel("czas lokalny")
    pl.ylabel("wysokość")
    pl.xticks(np.arange(min(czas), max(czas) + 2, 1))
    return fig

#wykres 2D wys(h)
czas = [i for i in range(24)] #posłuży jako oś X
fig = draw_2D(czas,h,"blue",'-')
draw_2D(czas,h1,"green",'--')
draw_2D(czas,h2,"red",':')




#wykres 3D xyz
fig2 = pl.figure("Wykres xyz Warszawa", figsize=(5,4))
ax2 = fig2.add_subplot(111,projection='3d')
ax2.scatter(x,y,z,c='blue',marker='o')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')


fig3 = pl.figure("Wykres xyz Dżakarta (fi = 0)", figsize=(5,4))
ax3 = fig3.add_subplot(111,projection='3d')
ax3.scatter(x1,y1,z1,c='green',marker='o')
ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.set_zlabel('Z')


fig4 = pl.figure("Wykres xyz Canberra", figsize=(5,4))
ax4 = fig4.add_subplot(111,projection='3d')
ax4.scatter(x2,y2,z2,c='red',marker='o')
ax4.set_xlabel('X')
ax4.set_ylabel('Y')
ax4.set_zlabel('Z')

ax3.ticklabel_format(axis="both", style="sci",scilimits = (0,0))
move_figure(fig, 700, 5)
move_figure(fig2, 100, 500)
move_figure(fig3, 700, 500)
move_figure(fig4, 1300, 500)
pl.show()
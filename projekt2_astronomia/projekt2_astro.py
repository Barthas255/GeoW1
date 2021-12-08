from astropy.time import Time
from astropy.coordinates import Angle
from astropy import units as u
import matplotlib.pyplot as pl
import matplotlib
import numpy as np

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
    return z, Az


#Zamiana na xyz
def x_y_z(z, Az, r = 1):
    z = np.deg2rad(z)
    Az = np.deg2rad(Az)
    x = r*np.sin(z)*np.cos(Az)
    y = r * np.sin(z) * np.sin(Az)
    z = r * np.cos(z)
    return x, y, z


#przesuwanie utworzonego wykresu do punktu x y
def move_figure(f, x, y):
    backend = matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
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


#Wykres wysokości od czasu
def draw_2D(czas, h, clr, ls_type):
    fig = pl.figure("Wysokość w czasie", figsize=(7,5))
    pl.plot(czas, h, color=clr, ls = ls_type)
    pl.grid(True)
    pl.xlabel("czas lokalny")
    pl.ylabel("wysokość")
    pl.xticks(np.arange(min(czas), max(czas) + 2, 1))
    pl.gca().spines['top'].set_visible(False)
    pl.gca().spines['right'].set_visible(False)
    return fig


#wykres 2D h(t)
czas = [i for i in range(24)] #posłuży jako oś X
fig = draw_2D(czas,h,"blue",'-')
draw_2D(czas,h1,"green",'--')
draw_2D(czas,h2,"red",':')


#wykres 3D xyz

#Płn półkula
fig2 = pl.figure("Wykres xyz Warszawa", figsize=(5,4))
ax2 = fig2.add_subplot(111,projection='3d')
ax2.scatter(x,y,z,c='blue',marker='o', zorder=2)
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')

#Równik
fig3 = pl.figure("Wykres xyz Dżakarta (fi = 0)", figsize=(5,4))
ax3 = fig3.add_subplot(111,projection='3d')
ax3.scatter(x1,y1,z1,c='green',marker='o', zorder=2)
ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.set_zlabel('Z')

#Płd półkula
fig4 = pl.figure("Wykres xyz Canberra", figsize=(5,4))
ax4 = fig4.add_subplot(111,projection='3d')
ax4.scatter(x2,y2,z2,c='red',marker='o', zorder=2)
ax4.set_xlabel('X')
ax4.set_ylabel('Y')
ax4.set_zlabel('Z')

#upewnienie sie, ze wykres na równiku jest w notacji naukowej
ax3.ticklabel_format(axis="both", style="sci",scilimits = (0,0))

#przesunięcie figur na ekranie
move_figure(fig, 700, 5)
move_figure(fig2, 100, 550)
move_figure(fig3, 700, 550)
move_figure(fig4, 1300, 550)

#sfery referencyjne
r = 0.95
u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
a = r*np.cos(u) * np.sin(v)
b = r*np.sin(u) * np.sin(v)
c = r*np.cos(v)
ax2.plot_surface(a, b, c, cmap = 'winter', alpha=0.5, zorder=1)
ax3.plot_surface(a, b, c, cmap = 'summer', alpha=0.5, zorder=1)
ax4.plot_surface(a, b, c, cmap = 'autumn', alpha=0.5, zorder=1)

pl.show()

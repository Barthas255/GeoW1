import numpy as np
from shapely.geometry.polygon import Polygon

#przyklad
fi1 = 51.70102972777778
lam1 = 18.175462347222222

A = (50.25, 20.75)
B = (50.0, 20.75)
C = (50.25, 21.25)
D = (50.0, 21.25)
S = (50.125, 21.0)
M = (50.12525870, 21.00063676) #S i M policzone z zadania 3
P = (A, B, C, D, S, M)

a = 6378137
b = 6356752.3141
e2 = 0.00669438002290
eprim2 = 0.00673949677548

a0 = 1 - e2/4 - 3*(e2**2)/64 - 5*(e2**3)/256
a2 = 3/8 * (e2 + (e2**2)/4 + 15*(e2**3)/128)
a4 = 15/256 * (e2**2 + 3*(e2**3)/4)
a6 = 35*(e2**3)/3072


def gauss_kruger(fi, lam, lam0):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    lam0 = np.deg2rad(lam0)
    sigma = a*(a0*fi - a2*np.sin(2*fi) + a4*np.sin(4*fi) - a6*np.sin(6*fi))
    t = np.tan(fi)
    eta2 = eprim2*np.cos(fi)
    l = lam - lam0
    N = a / (np.sqrt(1 - e2*(np.sin(fi)**2)))
    x = sigma + (l**2)/2 * N * np.sin(fi) * np.cos(fi) * \
        (1 + (l**2)/12*(np.cos(fi)**2)*(5-t**2+9*eta2+4*(eta2**2)) + (l**4)/360*(np.cos(fi)**4)*(61-58*(t**2)+t**4+270*eta2-330*eta2*(t**2)))

    y = l * N * np.cos(fi) * \
        (1 + (l**2)/6*(np.cos(fi)**2)*(1-t**2+eta2) + (l**4)/120*(np.cos(fi)**4)*(5-18*(t**2)+t**4+14*eta2-58*eta2*(t**2)))
    return x, y


def przedzialy_lamdy(lam):
    if lam < 16.5:
        return 5, 15
    elif lam >= 16.5 and lam < 19.5:
        return 6, 18
    elif lam >= 19.5 and lam < 22.5:
        return 7, 21
    elif lam >= 22.5:
        return 8, 24
    return 0, 0


def gauss_kruger_odwrotny(x, y, lam, uklad):
    B1_prev = 0
    B1_next = x / (a*a0)
    eps = np.deg2rad(0.00001 / 3600)
    while abs(B1_next - B1_prev) > eps:
        B1_prev = B1_next
        B1_next = x / (a*a0) + a2/a0 * np.sin(2*B1_prev) - a4/a0 * np.sin(4*B1_prev) + a6/a0 * np.sin(6*B1_prev)

    B = B1_next
    eta2 = eprim2 * np.cos(B)
    t = np.tan(B)
    N = a / (np.sqrt(1 - e2 * (np.sin(B) ** 2)))
    M = a*(1-e2) / np.sqrt((1 - e2 * np.sin(B)**2)**3)
    R = np.sqrt(M*N)
    mgk = 1 + (y**2)/(2*R**2) + (y**4)/(24*R**4)
    zgk = mgk - 1
    pgk = mgk**2
    zpgk = pgk - 1
    '''fi = fi - t/2 * ((y/(N*m0))**2 * (1+eta2) - 1/12 * (y/(N*m0))**4 * (5+3*t**2+6*eta2-6*eta2*t**2-3*eta2**2-9*t**2*eta2**2)
                      + 1/360 * (y/(N*m0))**6 * (61 + 90*t**2 + 45*t**4 + 107*eta2 - 162*t**2*eta2 - 45*t**4*eta2))
    l = 1/np.cos(fi) * (y/(N*m0) - 1/6 * (y/(N*m0))**3 * (1+2*t**2+eta2)
                        + 1/120 * (y/N*m0)**5 * (5+28*t**2+24*t**4+6*eta2+8*eta2*t**2))''' #inny wzór, mniej dokladny
    fi = B - (y**2*t)/(2 * M * N) * (1 - (y**2)/(12 * N**2) * (5+3*t**2+eta2-9*eta2*t**2-4*eta2**2)
                                      + (y**4)/(360*N**4) * (61 + 90*t**2 + 45*t**4))
    l = y/(N*np.cos(B)) * (1 - y**2/(6*N**2) * (1+2*t**2+eta2)
                              + y**4/(120*N**4) * (5+28*t**2+24*t**4+6*eta2+8*eta2*t**2))
    fi = np.rad2deg(fi)
    l = np.rad2deg(l)
    if uklad == "1992":
        l = float(l)
        l += 19
        m_ukladu = 0.9993*mgk
        z_ukladu = m_ukladu - 1
        p_ukladu = 0.9993**2 * pgk
    elif uklad == "2000":
        lam0 = przedzialy_lamdy(lam)
        l = float(l)
        l += lam0[1]
        m_ukladu = 0.999923 * mgk
        p_ukladu = 0.999923 ** 2 * pgk

    z_ukladu = m_ukladu - 1
    zp_ukladu = p_ukladu - 1
    return fi, l, mgk, zgk, m_ukladu, z_ukladu,    pgk, p_ukladu, zpgk, zp_ukladu


def dziewiec_dwa(fi, lam):
    m0 = 0.9993
    xgk, ygk = gauss_kruger(fi, lam, 19)
    x = m0 * xgk - 5300000
    y = m0 * ygk + 500000
    return x, y


def dziewiec_dwa2gauss(x, y):
    m0 = 0.9993
    xgk = (x + 5300000) / m0
    ygk = (y - 500000) / m0
    return xgk, ygk


def dwa_tysiace(fi, lam):
    m0 = 0.999923
    x = 0.0
    y = 0.0
    nr, pld_osiowy = przedzialy_lamdy(lam)
    xgk, ygk = gauss_kruger(fi, lam, pld_osiowy)
    x = m0 * xgk
    y = m0 * ygk + nr * 1000000 + 500000
    return x, y


def dwa_tysiace2gauss(x,y,lam):
    m0 = 0.999923
    xgk = x / m0
    nr, pld_osiowy = przedzialy_lamdy(lam)
    ygk = (y - nr*1000000 - 500000) / m0
    return xgk, ygk


#zamiana stopni na stopnie,minuty,sekundy do 5 miejsca po przecinku
def stopnie(omega):
    omega2 = (omega - int(omega)) * 60
    omega3 = (omega2 - int(omega2)) * 60
    omega = int(omega)
    omega2 = int(omega2)
    omega3 = round(omega3, 5)
    ds = u'\N{DEGREE SIGN}'
    f = str(omega) + ds +str(omega2) + "'"+str(omega3) +"'' "
    return f


def pole_elipsoidalne(fi1, lam1, fi2, lam2):
    lam1 = np.deg2rad(lam1)
    lam2 = np.deg2rad(lam2)

    fi1 = np.deg2rad(fi1)
    fi2 = np.deg2rad(fi2)
    e = e2**0.5
    return  abs(
            (b**2 * (lam2 - lam1) / 2) * (((np.sin(fi2) /
            (1 - e2 * (np.sin(fi2)**2))) + (1 / (2 * e))
            * np.log(
                    (1 + e * np.sin(fi2)) / (1 - e * np.sin(fi2))))
                    - ((np.sin(fi1) / (1 - e2 * (np.sin(fi1)**2))) + (1 / (2 * e)) *
                np.log(
                    (1 + e * np.sin(fi1)) / (1 - e * np.sin(fi1)))))
            )


#test
'''print("Dane: " ,fi1,lam1)
print("gauss dla ",gauss_kruger(fi1,lam1,19))
x92,y92 = dziewiec_dwa(fi1,lam1)
x00,y00 = dwa_tysiace(fi1,lam1)
print("1992: ", x92,y92)
print("2000: ", x00,y00)
x92,y92 = dziewiec_dwa2gauss(x92,y92)
x00,y00 = dwa_tysiace2gauss(x00,y00,lam1)
print("1992 > Gauss: ", x92,y92)
print("2000 > Gauss: ", x00,y00)
print("1992 > Gauss > Geo: ", gauss_kruger_odwrotny(x92,y92,lam1,"1992"))
print("2000 > Gauss > Geo: ", gauss_kruger_odwrotny(x00,y00,lam1,"2000"))'''

slownik = {1:"A",2:"B",3:"C",4:"D",5:"S",6:"M"}
pola = []
pelipsy = []
pgk = []
p2000 = []
p1992 = []
pola.append((pelipsy,pgk,p2000,p1992))
pom = 0
for i in P:
    pom += 1
    print("Punkt:",slownik[pom], stopnie(i[0]), stopnie(i[1]))
    xgk, ygk = gauss_kruger(i[0], i[1], 19)
    print("Xgk, Ygk: ", round(xgk, 3), round(ygk, 3))
    x92, y92 = dziewiec_dwa(i[0], i[1])
    x00, y00 = dwa_tysiace(i[0], i[1])
    if pom <= 4:
        pelipsy.append((i[0],i[1]))
        pgk.append((xgk, ygk))
        p1992.append((x92, y92))
        p2000.append((x00, y00))
    print("1992: ", round(x92, 3), round(y92, 3))
    print("2000: ", round(x00, 3), round(y00, 3))
    x92, y92 = dziewiec_dwa2gauss(x92, y92)
    x00, y00 = dwa_tysiace2gauss(x00, y00, i[1])
    print("1992 > Gauss: ", round(x92, 3), round(y92, 3))
    print("2000 > Gauss: ", round(x00, 3), round(y00, 3))
    x92, y92, m, z, m92, z92, pgk2, p92, zgk, zp92 = gauss_kruger_odwrotny(x92, y92, i[1], "1992")
    x00, y00, foo, bar, m00, z00, nothing, p00, something, zp00 = gauss_kruger_odwrotny(x00, y00, i[1], "2000")
    print("1992 > Gauss > Geo: ", x92, y92)
    print("2000 > Gauss > Geo: ", x00, y00)
    print("mgk: ", round(m, 6))
    print("Kgk: ", round(z*1000, 2), "cm/1km")
    print("m92: ", round(m92, 6))
    print("z92: ", round(z92 * 1000, 2), "cm/1km")
    print("m00: ", round(m00, 6))
    print("z00: ", round(z00 * 1000, 2), "cm/1km")
    print('\n')
    print("pgk: ", round(pgk2, 6))
    print("Kgk(1ha): ", round(zgk*10000, 2), "cm/1ha")
    print("p92: ", round(p92, 6))
    print("z92: ", round(zp92 * 10000, 2), "cm/1ha")
    print("p00: ", round(p00, 6))
    print("z00: ", round(zp00 * 10000, 2), "cm/1ha")
    print('\n --------------------------------------------- \n')
print("Pole elipsoidalne: ", round(pole_elipsoidalne(pola[0][0][0][0], pola[0][0][0][1], pola[0][0][3][0], pola[0][0][3][1])/1000000, 12), "km^2")
pola2 = [(pola[0][1][0],pola[0][1][1],pola[0][1][3],pola[0][1][2]), (pola[0][2][0],pola[0][2][1],pola[0][2][3],pola[0][2][2]),
         (pola[0][3][0],pola[0][3][1],pola[0][3][3],pola[0][3][2])]
slownik = {1:"Pole GK:",2:"Pole 1992:",3:"Pole 2000:"}
pom = 1
for i in pola2:
    geom = Polygon(i)
    print(slownik[pom],round(geom.area/1000000, 12),"km^2")
    pom += 1
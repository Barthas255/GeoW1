import numpy as np

#GRS80
e2 = 0.00669437999013
a = 6378137
H = 100 #zalozenie ze h = 100

#Krasowski
akr = 6378245
bkr = 6356863
e2kr = (akr**2 - bkr**2) / (akr**2)

def geo2xyz(fi, lam, h): #GRS 80
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a/np.sqrt(1-e2*np.sin(fi)**2)

    x = (N+h)*np.cos(fi)*np.cos(lam)
    y = (N+h)*np.cos(fi)*np.sin(lam)
    z = (N*(1-e2)+h)*np.sin(fi)

    return x, y, z

def hirvonen(x,y,z):
    #1.
    r = np.sqrt(x**2 + y**2)
    #2.
    fi_next = np.arctan(z/r * (1/(1-e2kr)))
    fi_prev = 0
    #3-5.
    while abs(fi_next - fi_prev) > np.deg2rad(0.00005 / 3600):
        fi_prev = fi_next
        N = akr / np.sqrt(1 - e2kr * (np.sin(fi_next) ** 2))
        h = r / np.cos(fi_next) - N
        fi_next = np.arctan(z/r / (1-e2kr * (N/(N+h))))
    #6
    lam = np.arctan(y/x)
    N = akr / np.sqrt(1 - e2kr * np.sin(fi_next) ** 2)
    h = r/np.cos(fi_next) - N
    #print("współrzędne przed: ", x, y, z)
    x = (N+h)*np.cos(fi_next)*np.cos(lam)
    y = (N+h)*np.cos(fi_next)*np.sin(lam)
    z = (N*(1-e2kr)+h)*np.sin(fi_next)
    #print("współrzędne po: ", x, y, z)
    lam = np.rad2deg(lam)
    fi_next = np.rad2deg(fi_next)

    return fi_next,lam,h


def transformacja_bursy_wolfa(xp,yp,zp):
    x0 = -33.4297
    y0 = 146.5746
    z0 = 76.2865
    alfa = -0.35867 #sekundy
    beta = -0.05283 #sekundy
    gamma = 0.84354 #sekundy
    kappa = 0.0000008407728

    alfa = np.deg2rad(alfa / 3600)
    beta = np.deg2rad(beta / 3600)
    gamma = np.deg2rad(gamma / 3600)

    xyzp = np.array([[xp],[yp],[zp]])
    abgk = np.array([kappa,gamma,-beta, -gamma,kappa,alfa, beta,-alfa,kappa])
    abgk = abgk.reshape((3,3))
    xyz0 = np.array([[x0],[y0],[z0]])
    xyzw = xyzp + (abgk @ xyzp) + xyz0
    return xyzw[0][0],xyzw[1][0],xyzw[2][0]


#zamiana tekstowa stopni
def stopnie(omega):
    omega2 = (omega - int(omega)) * 60
    omega3 = (omega2 - int(omega2)) * 60
    omega = int(omega)
    omega2 = int(omega2)
    omega3 = round(omega3, 5)
    ds = u'\N{DEGREE SIGN}'
    f = str(omega) + ds +str(omega2) + "'"+str(omega3) +"'' "
    return f


#dzialajacy przyklad
'''fii = 50.00
lamd = 20.00
H1 = 100
fii_text = stopnie(fii)
lamd_text = stopnie(lamd)
print("fi,lambda,h  grs80: " ,fii_text,lamd_text,H1)
raz,dwa,trzy = geo2xyz(fii,lamd,H1)
raz1,dwa2,trzy3= transformacja_bursy_wolfa(raz,dwa,trzy)
print("transformacja bursy wolfa: ",round(raz1,3),round(dwa2,3),round(trzy3,3))
raz1,dwa2,trzy3 = hirvonen(raz1,dwa2,trzy3)
raz1 = stopnie(raz1)
dwa2 = stopnie(dwa2)
print("hirvonen: ",raz1,dwa2,round(trzy3,3))'''

#punkty z poprzedniego zadania
A = (50.25, 20.75)
B = (50.0, 20.75)
C = (50.25, 21.25)
D = (50.0, 21.25)
S = (50.125, 21.0)
M = (50.12525870, 21.00063676) #S i M policzone z zadania 3
P = (A, B, C, D, S, M)
slownik = {1:"A", 2:"B", 3:"C", 4:"D", 5:"S", 6:"M"}
tab = []
tab2 = []
for i in range(6):
    print("Współrzędne punktu " + slownik[i+1] + " na elipsoidzie GRS80: ", stopnie(P[i][0]), stopnie(P[i][1]), H, "m")
print("-" * 85)
for i in range(6):
    x, y, z = geo2xyz(P[i][0], P[i][1], H)
    tab.append((x,y,z))
    print("Współrzędne punktu " + slownik[i+1] + " na elipsoidzie GRS80: ", round(x, 3), round(y, 3), round(z, 3))
print("-"*85)
for i in range(6):
    xkr, ykr, zkr = transformacja_bursy_wolfa(*tab[i])
    tab2.append((xkr,ykr,zkr))
    print("Współrzędne punktu " + slownik[i+1] + " na elipsoidzie Krasowskiego: ", round(xkr, 3), round(ykr, 3), round(zkr, 3))
print("-"*85)
for i in range(6):
    phi, lam, h = hirvonen(*tab2[i])
    phi = stopnie(phi)
    lam = stopnie(lam)
    print("Współrzędne geodezyjne punktu " + slownik[i+1], phi, lam, round(h, 3),"m")
print("-"*85)
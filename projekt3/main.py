import numpy as np

nr = 0

fiA=50.25 + nr * 0.25
lambdaA=20.75
fiC=50.25 + nr * 0.25
lambdaC=12.25
fiB=50.0 + nr * 0.25
lambdaB=20.75
fiD=50.0 + nr * 0.25
lambdaD=21.25


a = 6378137  # metry
e2 = 0.0066943800290  # bez jednostek
b = a * (1 - e2) ** 0.5

def vincenty(fi_a, lambda_a, fi_b, lambda_b, a, e2):
    lambda_a = np.deg2rad(lambda_a)
    fi_a = np.deg2rad(fi_a)
    lambda_b = np.deg2rad(lambda_b)
    fi_b = np.deg2rad(fi_b)

    #1.
    b = a * np.sqrt(1 - e2)
    f = 1 - b/a
    #2.
    deltaLambda = lambda_b - lambda_a
    #3.
    Ua = np.arctan((1 - f) * np.tan(fi_a))
    Ub = np.arctan((1 - f) * np.tan(fi_b))
    #4.
    L_before = 0
    L_after = deltaLambda
    cos2alfa, sinDelta, cosDelta, cos_podwojny_delta_m, delta = 0,0,0,0,0
    #5-12
    while np.abs(L_after - L_before) > 0.0000000000027:
        #5.
        L_before = L_after
        L = L_before
        sinDelta = np.sqrt((np.cos(Ub)*np.sin(L))**2 + (np.cos(Ua)*np.sin(Ub) - np.sin(Ua)*np.cos(Ub)*np.cos(L))**2)
        #6.
        cosDelta = np.sin(Ua)*np.sin(Ub) + np.cos(Ua)*np.cos(Ub)*np.cos(L)
        #7.
        delta = np.arctan(sinDelta/cosDelta)
        #8.
        sin_alfa = np.cos(Ua)*np.cos(Ub)*np.sin(L)/sinDelta
        #9.
        cos2alfa = 1 - sin_alfa**2
        #10.
        cos_podwojny_delta_m = cosDelta - 2*np.sin(Ua)*np.sin(Ub)/cos2alfa
        #11.
        C = f/16 * cos2alfa * (4 + f*(4 - 3*cos2alfa))
        #12.
        L_after = deltaLambda + (1 - C) * f * sin_alfa * (delta + C * sinDelta * (cos_podwojny_delta_m + C*cosDelta * (-1 + 2 * (cos_podwojny_delta_m**2))))

    L = L_after
    #13.
    u2 = (a**2 - b**2) * cos2alfa / (b**2)
    #14.
    A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175*u2)))
    #15.
    B = (u2/1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    #16.
    delta_delta = B * sinDelta * \
                  (cos_podwojny_delta_m + 1/4 * B *
                   (cosDelta * (-1 + 2*(cos_podwojny_delta_m**2)) - 1/6 * B * cos_podwojny_delta_m * (-3 + 4*(sinDelta**2)) * (-3 + 4*(cos_podwojny_delta_m**2))))
    #17.
    Sab = b * A * (delta - delta_delta)
    #18-19
    licznik_ab = np.cos(Ub) * np.sin(L)
    mianownik_ab = np.cos(Ua) * np.sin(Ub) - np.sin(Ua) * np.cos(Ub) * np.cos(L)
    licznik_ba = np.cos(Ua) * np.sin(L)
    mianownik_ba = -np.sin(Ua) * np.cos(Ub) + np.cos(Ua) * np.sin(Ub) * np.cos(L)
    Aab = np.arctan(licznik_ab/mianownik_ab)
    Aba = np.arctan(licznik_ba/mianownik_ba) + np.pi
    Aab = np.rad2deg(Aab)
    Aba = np.rad2deg(Aba)

    if mianownik_ab < 0: Aab += 180
    elif licznik_ab < 0: Aab += 360
    if Aab > 360:
        Aab -= 360

    if mianownik_ba < 0:
        Aba += 180
    elif licznik_ba < 0:
        Aba += 360

    if Aba > 360:
        Aba -= 360
    if Aab > 360:
        Aab -= 360
    return np.abs(Sab), Aab, Aba


def kivioj(a,e2,prev_fi, prev_lamda, ds, prev_Az12):
    prev_fi = np.deg2rad(prev_fi)
    prev_lamda = np.deg2rad(prev_lamda)
    prev_Az12 = np.deg2rad(prev_Az12)

    M = (a*(1-e2)) / (np.sqrt(((1-e2*(np.sin(prev_fi))**2)**3))) #glowne promienie krzywizny
    N = a/(np.sqrt(1-e2*(np.sin(prev_fi))**2)) #glowne promienie krzywizny
    delta_fi = (ds*np.cos(prev_Az12))/M #przyblizenie przyrostu szerokosci
    fi_sr_ds = prev_fi + 1/2 * delta_fi #srednia wartosc szerokosci elementu ds

    A12_sr = prev_Az12 + ds/N * np.sin(prev_Az12)*np.tan(fi_sr_ds) #sredni azymut
    #lepsze przyblizenia przyrostu wartosci fi, lambda, azymut
    delta_fi = ds*np.cos(A12_sr)/M
    delta_lambda = (ds*np.sin(A12_sr))/(N*np.cos(fi_sr_ds))
    delta_Az12 = (ds * np.sin(A12_sr)*np.tan(fi_sr_ds))/N

    next_fi = prev_fi + delta_fi
    next_lamda = prev_lamda + delta_lambda
    next_Az12 = prev_Az12 + delta_Az12

    next_fi = np.rad2deg(next_fi)
    next_lamda = np.rad2deg(next_lamda)
    next_Az12 = np.rad2deg(next_Az12)

    return next_fi,next_lamda,next_Az12


def pole(fi1, lam1, fi2, lam2, e2, b):
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


s,az12,az21 = vincenty(fiA,lambdaA,fiD,lambdaD,a,e2)
print("dlugosc odcinka AD: " +str(s) +" m")

fi = fiA
lam = lambdaA
ds = 1000
s = s/2
while s > 1000:
    fi,lam,az12 = kivioj(a,e2,fi,lam,ds,az12)
    s -= 1000

if s > 0:
    ds = s
    fi,lam,az12 = kivioj(a,e2,fi,lam,ds,az12)

fi_sr_szer = (fiA+fiD)/2
lam_sr_szer = (lambdaA+lambdaD)/2

odl,Az12,Az21 = vincenty(fi,lam,fi_sr_szer,lam_sr_szer,a,e2)
print("punkt średniej szerokości: " + str(fi_sr_szer) +" " + str(lam_sr_szer), "\npunkt środkowy: " +str(fi) +" "+str(lam),
      "\nróżnica odległości punktu środkowego i średniej szerokości: " + str(odl),
      "\nazymut wprost pkt środkowego: " +str(Az12), "\nazymut odwrotny pkt środkowego: " +str(Az21) ,
      "\npole prostokąta: " + str(pole(fiA,lambdaA,fiD,lambdaD,e2,b)) + " m^2")

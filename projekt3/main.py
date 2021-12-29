import numpy as np
'''
Oznaczenia

a,b - dluzsza i krotsza półoś elipsoidy
e, f - pierwszy mimośród, spłaszczenie elipsoidy
fi,lambda - szer. i dł. geodezyjna
deltalambda - roznica dlugosci geodezyjnej
Aab, Aba - azymur prosty i odwrotny
alfa - azymut linii geodezyjnej na równiku
U - szerokość zredukowana
L - różnica długości na sferze pomocniczej
delta - odległość kątowa pomiędzy punktami na sferze
deltam - odległość kątowa na sferze od równika do punktu środkowego linii geodezyjnej

Dane
fiA, lambdaA - wsp geodezyjne punktu A
fiB, lambdaB - wsp geodezyjne punktu B
a, e^2 - parametry elipsoidy

szukane
Sab - dlugosc linii geodezyjnej pomiędzy punktami A i B
Aab, Aba - azymut prosty i odwrotny linii geodezyjnej
'''

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

    if licznik_ba > 0 and mianownik_ba > 0:
        Aba += 180
    elif mianownik_ba < 0:
        Aba += 360
    elif licznik_ab < 0:
        Aba += 540

    if Aba > 360:
        Aba -= 360
    if Aab > 360:
        Aab -= 360
    return np.abs(Sab), Aab, Aba


def kivioj(a,e2,prev_fi, prev_lamda, ds, prev_Az12):
    M = (a*(1-e2)) / (np.sqrt(((1-e2*(np.sin(prev_fi))**2)**3))) #glowne promienie krzywizny
    N = a/(np.sqrt(1-e2*(np.sin(prev_fi))**2)) #glowne promienie krzywizny
    delta_fi = (ds*np.cos(prev_Az12))/M #przyblizenie przyrostu szerokosci
    fi_sr_ds = prev_fi + 1/2 * delta_fi #srednia wartosc szerokosci elementu ds

    sr_M = (a * (1 - e2)) / (np.sqrt(((1 - e2 * (np.sin(fi_sr_ds)) ** 2) ** 3))) #srednie M
    sr_N = a / (np.sqrt(1 - e2 * (np.sin(fi_sr_ds)) ** 2)) #srednie N
    A12_sr = prev_Az12 + ds/sr_N * np.sin(prev_Az12)*np.tan(fi_sr_ds) #sredni azymut
    #lepsze przyblizenia przyrostu wartosci fi, lambda, azymut
    delta_fi = ds*np.cos(A12_sr)/sr_M
    delta_lambda = (ds*np.sin(delta_fi))/(sr_N*np.cos(fi_sr_ds))
    delta_Az12 = ds/sr_N * np.sin(A12_sr)*np.tan(fi_sr_ds)

    next_fi = prev_fi + delta_fi
    next_lamda = prev_lamda + delta_lambda
    next_Az12 = prev_Az12 + delta_Az12

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
    az21 = az12 + 180
    s -= 1000

if s > 0:
    ds = s
    fi,lam,az12 = kivioj(a,e2,fi,lam,ds,az12)
    az21 = az12 + 180

fi_sr_szer = (fiA+fiD)/2
lam_sr_szer = (lambdaA+lambdaD)/2

odl,foo,bar = vincenty(fi,lam,fi_sr_szer,lam_sr_szer,a,e2)
print("punkt średniej szerokości: " + str(fi_sr_szer) +" " + str(lam_sr_szer), "\npunkt środkowy: " +str(fi) +" "+str(lam),
      "\nróżnica odległości punktu środkowego i średniej szerokości: " + str(odl),
      "\nazymut wprost pkt środkowego: " +str(az12), "\nazymut odwrotny pkt środkowego: " +str(az21) ,
      "\npole prostokąta: " + str(pole(fiA,lambdaA,fiD,lambdaD,e2,b)) + " m^2")










'''
from math import *

# parametry dla grs 80
a = 6378137
e2 = 0.00669437999013
b = a * (1 - e2)**0.5
f = 1 / 298.257222101

def M(phi):
    return (a*(1 - e2)) / ((1 - e2*sin(radians(phi))**2)**3)**0.5

def N(phi):
    return a / (1 - e2*sin(radians(phi))**2)**0.5

def pole(lam1, phi1, lam2, phi2):
    lam1 = radians(lam1)
    lam2 = radians(lam2)

    phi1 = radians(phi1)
    phi2 = radians(phi2)
    e = e2**0.5
    return  abs(
            (b**2 * (lam2 - lam1) / 2) * (((sin(phi2) /
            (1 - e2 * (sin(phi2)**2))) + (1 / (2 * e))
            * log(
                    (1 + e * sin(phi2)) / (1 - e * sin(phi2))))
                    - ((sin(phi1) / (1 - e2 * (sin(phi1)**2))) + (1 / (2 * e)) *
                log(
                    (1 + e * sin(phi1)) / (1 - e * sin(phi1)))))
            )

class Kivioj():
    def __init__(self, S, Az_p, phi_p, lam_p):
        ds = 1250
        n = int(S / ds)
        ds_o = S % ds
        for i in range(n):
            d_phi = (ds * cos(Az_p)) / M(phi_p)
            d_Az = (sin(Az_p) * tan(radians(phi_p)) * ds) / N(phi_p)

            Sr_phi = phi_p + degrees(d_phi) / 2
            Sr_az = Az_p + d_Az / 2
            d_phi_2 = (ds * cos(Sr_az)) / M(Sr_phi)
            d_lam = (ds * sin(Sr_az)) / (N(Sr_phi) * cos(radians(Sr_phi)))

            d_Az2 = (sin(Sr_az) * tan(radians(Sr_phi)) * ds) / N(Sr_phi)

            phi_p = phi_p + degrees(d_phi_2)
            lam_p = lam_p + degrees(d_lam)
            Az_p = Az_p + d_Az2

        d_phi = (ds_o * cos(Az_p)) / M(phi_p)
        d_Az = (sin(Az_p) * tan(radians(phi_p)) * ds_o) / N(phi_p)
        Sr_phi = phi_p + degrees(d_phi / 2)

        Sr_az = Az_p + d_Az / 2
        d_phi_2 = (ds_o * cos(Sr_az)) / M(Sr_phi)
        d_lam = (ds_o * sin(Sr_az)) / (N(Sr_phi) * cos(radians(Sr_phi)))

        d_Az2 = (sin(Sr_az) * tan(radians(Sr_phi)) * ds_o) / N(Sr_phi)
        phi_k = phi_p + degrees(d_phi_2)
        lam_k = lam_p + degrees(d_lam)

        Az = Az_p + d_Az2

        self.result = {
            "phi": phi_k,
            "lambda": lam_k,
            "Az": Az
        }


class Vincent():

    def __init__(self, lambda1, lambda2, fi_1, fi_2):
        U1 = atan((1 - f) * tan(radians(fi_1)))
        U2 = atan((1 - f) * tan(radians(fi_2)))

        L = radians(lambda2 - lambda1)

        difference = 100
        precision = radians(0.000001 / 3600)

        while difference >= precision:
            sins = ((cos(U2) * sin(L)) ** 2 + (
                    cos(U1) * sin(U2) - sin(U1) *
                    cos(U2) * cos(L)) ** 2) ** 0.5

            coss = sin(U1) * sin(U2) + cos(U1) * cos(U2) * cos(L)
            sig = atan(sins / coss)

            sina = (cos(U1) * cos(U2) * sin(L)) / sins
            cos2a = 1 - sina ** 2

            cos2sm = coss - ((2 * sin(U1) * sin(U2)) / cos2a)
            C = (f / 16) * cos2a * (4 + f * (4 - 3 * cos2a))

            prev_lambda = radians(lambda2 - lambda1) + (1 - C) * f * sina * (
                    sig + C * sins *
                    (cos2sm + C * coss * ((-1) +
                                          2 * cos2sm ** 2)))

            difference = abs(prev_lambda - L)
            L = prev_lambda

        u2 = ((a ** 2 - b ** 2) / b ** 2) * cos2a

        A = 1 + (u2 / 16384) * (4096 + u2 *
                                ((-768) + u2 * (320 - 175 * u2)))

        B = (u2 / 1024) * (256 + u2 *
                           ((-128) + u2 * (74 - 47 * u2)))

        d_sig = B * sins * (cos2sm + (1 / 4) * B * (
                coss * ((-1) + 2 * cos2sm ** 2) - (1 / 6) * B * cos2sm * ((-3) + 4 * sins ** 2) * (
                (-3) + 4 * cos2sm ** 2)))

        # odległość
        s = b * A * (sig - d_sig)
        self.s = s

        # azymuty
        y = cos(U2) * sin(prev_lambda)
        x = cos(U1) * sin(U2) - sin(U1) * cos(U2) * cos(prev_lambda)


        # poprawa azymutów
        if (y > 0 and x > 0):
            Az12 = atan(y / x)
        elif (y > 0 and x < 0):
            Az12 = atan(y / x) + pi
        elif (y < 0 and x < 0):
            Az12 = atan(y / x) + pi
        elif (y < 0 and x > 0):
            Az12 = atan(y / x) + 2 * pi

        self.Az12 = Az12
        y = cos(U1) * sin(prev_lambda)
        x = (-sin(U1)) * cos(U2) + cos(U1) * sin(U2) * cos(prev_lambda)
        # poprawa azymutó dla odwrotnego
        if (y > 0 and x > 0):
            Az21 = atan(y / x) + pi
        elif (y > 0 and x < 0):
            Az21 = atan(y / x) + 2 * pi
        elif (y < 0 and x < 0):
            Az21 = atan(y / x) + 2 * pi
        elif (y < 0 and x > 0):
            Az21 = atan(y / x) + 3 * pi

        self.Az21 = Az21


if __name__ == "__main__":
    a_fi = 50.25
    a_lam = 20.75
    b_fi = 50.25
    b_lam = 21.25
    c_fi = 50.
    c_lam = 20.75
    d_fi = 50.
    d_lam = 21.25



    a_d = Vincent(a_lam,d_lam, a_fi, d_fi)
    print("odleglosc AD i Azymut wprost: ",a_d.s, 'm', ',',  degrees(a_d.Az12), "°")
    middle_fi = (a_fi + d_fi) /2
    middle_lambda = (a_lam + d_lam) /2


    print("Punkt średniej szerokości: " +
          str((middle_fi)) + "°"
          + str((middle_lambda)) + "°"
          )
    # punkt środkowy

    kiv = Kivioj(a_d.s/2, a_d.Az12, a_fi, a_lam)
    print("Punkt środkowy AD: "
          + str((kiv.result['phi'])) + "°"
          + str((kiv.result['lambda'])) + "°"

          )
    # pole powieszchni
    print("Pole powierzchni czworokata: "
          + str(pole(a_lam, a_fi, d_lam, d_fi)), 'm^2')


    srod_fi = kiv.result['phi']
    srod_lam = kiv.result['lambda']


    roznica = Vincent(middle_lambda,srod_lam, middle_fi, srod_fi)
    print("różnica odległości między punktem środkowym a średniej szerokości", roznica.s, 'm')
    print("Azymut wprost i odwrotny punktu średniej szerokości i punktu środkowego",
          degrees(roznica.Az12) , "° ", degrees(roznica.Az21), "°")


    #def __init__(self, lambda1, lambda2, fi_1, fi_2):
    # test = Vincent(20.75,20.75, 50.25, 50)
    # print(test.s)'''




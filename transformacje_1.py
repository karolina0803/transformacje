# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 17:23:08 2024

@author: USER
"""
import math as m

class Transformacje:
    
    def __init__(self, model: str = 'wgs84'):
        """
        Parametry
        ----------
        model : str, model elipsoidy, domyslny model to WGS84,
            dostepne modele:
                + WGS84 -> 'wgs84'
                + GRS80 -> 'grs80'
                + ELipsoida Krasowskiego -> 'elipsoida_krasowskiego'
        ----------
        Parametry elipsoid:
            a - duza polos elipsody - promien rownikowy
            b - mala polos elipsoidy - promien poludnikowy
            flat - splaszczenie elipsoidy
            e2 - kwadrat mimosrodu elipsoidy
        """
        if model == "wgs84":
            self.a = 6378137.0
            self.b = 6356752.31424518
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "elipsoida_krasowskiego":
            self.a = 6378245.0
            self.b = 6356863.019
        else:
            raise NotImplementedError(f"model {model} nie jest obslugiwany")
            
        self.flat = (self.a - self.b) / self.a
        self.e = m.sqrt(2 * self.flat - self.flat ** 2)
        self.e2 = (2 * self.flat - self.flat ** 2) 

    
    def obliczenie_N(self, phi):
        """
        funkcja pomocnicza - obliczenie promienia krzywizny w I wertykale
        """
        N = self.a/((1 - self.e2 * m.sin(phi)**2)**(1/2))
        return N
    
    def dlugosc_luku_poludnika(self, phi):
        '''
        funkcja obliczajaca dlugosc luku poludnika na podstawie phi punktu
        '''
        A0 = 1 - (self.e2/4) - ((3*(self.e2**2))/64) - ((5*(self.e2**3))/256)
        A2 = (3/8) * (self.e2 + ((self.e2**2)/4) + ((15*(self.e2**3))/128))
        A4 = (15/256) * (self.e2**2 + ((3*(self.e2**3))/4))
        A6 = (35*(self.e2**3))/3072
        
        # dlugosc luku poludnika
        sigma = self.a * (A0*phi - A2*m.sin(2*phi) + A4*m.sin(4*phi) - A6*m.sin(6*phi))
        
        return sigma
    
    def rad_degrees(self, val_rad):
        decimal_degrees = val_rad * (180 / m.pi)
        return decimal_degrees

    def rad_dms(self, val_rad):
        decimal_degrees = val_rad * (180 / m.pi)
        deg = m.trunc(decimal_degrees)
        dec_mins = (decimal_degrees - deg) * 60
        mins = m.trunc(dec_mins)
        sec = (dec_mins - mins) * 60
        return f'{deg:.0f}\xb0{mins:.0f}\'{sec:.5f}\"'

    def dms(self, val_dec_degrees):
        deg = m.trunc(val_dec_degrees)
        dec_mins = (val_dec_degrees - deg) * 60
        mins = m.trunc(dec_mins)
        sec = (dec_mins - mins) * 60
        return f'{deg:.0f}\xb0{mins:.0f}\'{sec:.5f}\"'

    def XYZ_philamh(self, X, Y, Z, output = 'dec'):
        '''
        funkcja wykorzystujaca algorytm hirvonena do przeliczenia wspolrzednych
        ortokartezjanskich (X,Y,Z) na wspolrzedne geodezyjne (dlugosc, szerokosc i wysokosc elipsoidalna)
        w wyniku iteracji otrzymujemy wspolrzedne z dokladnoscia ok 1 cm.
        
        -----
        parametry:
            X, Y, Z [float] - wspolrzedne w ukladzie ortokartezjanskim,
            
        -----
        returns:
            phi [str] - szerokosc geodezyjna
            lam [str] - dlugosc geodezyjna
            h [str] - wysokosc elipsoidalna [m]
            
        w zaleznosci od parametru output szerokosc i dlugosc geodezyjna
        moga byc podane w dwoch formatach:
                dec - stopnie dziesietne (domyslne)
                dms - stopnie, minuty, sekundy
                
        '''
        p = m.sqrt(X**2 + Y**2)
        phi = m.atan(Z/(p*(1-self.e2)))
        while True:
            phi_poprzednie = phi
            # N = self.a/((1 - self.e2 * m.sin(phi)**2)**(1/2))
            N = self.obliczenie_N(phi_poprzednie)
            h = p/m.cos(phi_poprzednie) - N
            phi = m.atan(Z/p * ((N*(1-self.e2)+h)/(N+h))**(-1))
            if abs(phi_poprzednie - phi)<(0.000001/206265):
                break
        N = self.obliczenie_N(phi)
        # N = self.a/((1 - self.e2 * m.sin(phi)**2)**(1/2))
        # h = p/m.cos(phi_poprzednie) - N
        h = p/m.cos(phi) - N
        lam = m.atan2(Y,X)
        
        if output == "dec":
             return f"{self.rad_degrees(phi)}, {self.rad_degrees(lam)}, {h:.3f}m"
        elif output == "dms":
             return f"{self.rad_dms(phi)}, {self.rad_dms(lam)}, {h:.3f}m"
    


            
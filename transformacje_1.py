# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 17:23:08 2024

@author: USER
"""
import math as m
import numpy as np
import sys 

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
        funkcja pomocnicza - obliczaja dlugosc luku poludnika na podstawie 
        szerokosci geodezyjnej punktu
        '''
        A0 = 1 - (self.e2/4) - ((3*(self.e2**2))/64) - ((5*(self.e2**3))/256)
        A2 = (3/8) * (self.e2 + ((self.e2**2)/4) + ((15*(self.e2**3))/128))
        A4 = (15/256) * (self.e2**2 + ((3*(self.e2**3))/4))
        A6 = (35*(self.e2**3))/3072
        
        # dlugosc luku poludnika
        sigma = self.a * (A0*phi - A2*m.sin(2*phi) + A4*m.sin(4*phi) - A6*m.sin(6*phi))
        
        return sigma
    
    def rad_degrees(self, val_rad):
        '''
        funkcja zamieniajaca wartosc w radianach na stopnie dziesietne
        '''
        decimal_degrees = val_rad * (180 / m.pi)
        return decimal_degrees

    def rad_dms(self, val_rad):
        '''
        funkcja zamieniajaca wartosc w radianach na stopnie minuty i sekundy 
        (jako lancuch znakow)
        '''
        decimal_degrees = val_rad * (180 / m.pi)
        deg = m.trunc(decimal_degrees)
        dec_mins = (decimal_degrees - deg) * 60
        mins = m.trunc(dec_mins)
        sec = (dec_mins - mins) * 60
        return f'{deg:.0f}\xb0{mins:.0f}\'{sec:.5f}\"'

    def dms(self, val_dec_degrees):
        '''
        funkcja zamieniajaca stopnie dziesietne na stopnie minuty i sekundy
        (jako lancuch znakow)
        '''
        deg = m.trunc(val_dec_degrees)
        dec_mins = (val_dec_degrees - deg) * 60
        mins = m.trunc(dec_mins)
        sec = (dec_mins - mins) * 60
        return f'{deg:.0f}\xb0{mins:.0f}\'{sec:.5f}\"'
    
    def dms_rad(self, val_dms):
        '''
        funkcja zamieniajaca wartosc stopnie minuty sekundy na radiany
        '''
        deg, mins, sec = val_dms
        rad = (deg + mins/60 + sec/3600)* (m.pi / 180)
        return rad
    
    def degrees_rad(self, val_dec_degrees):
        '''
        funkcja zamieniajaca wartosc stopni dziesietnych na radiany
        '''
        rad = val_dec_degrees * (m.pi / 180)
        return rad


    def XYZ_philamh(self, X, Y, Z, output = 'dec'):
        '''
        funkcja wykorzystujaca algorytm hirvonena do przeliczenia wspolrzednych
        ortokartezjanskich (X,Y,Z) na wspolrzedne geodezyjne (dlugosc, szerokosc i wysokosc elipsoidalna)
        w wyniku iteracji otrzymujemy wspolrzedne z dokladnoscia ok 1 cm.
        
        -----
        parametry:
            X, Y, Z [float] - wspolrzedne w ukladzie ortokartezjanskim
            
        -----
        returns:
            phi [float] - szerokosc geodezyjna
            lam [float] - dlugosc geodezyjna
            h [float] - wysokosc elipsoidalna [m]
            
        w zaleznosci od parametru output szerokosc i dlugosc geodezyjna
        moga byc podane w dwoch formatach:
                dec - stopnie dziesietne (wartosc domyslna)
                dms - stopnie, minuty, sekundy
                calc - zwrocone wartosci w radianach, np. jesli nie potrzeba
                    dokladnie tych wartosci, tylko do dalszych obliczen
                
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
        h = p/m.cos(phi) - N
        lam = m.atan2(Y,X)
        
        if output == "dec":
            return self.rad_degrees(phi), self.rad_degrees(lam), h
             
        elif output == "dms":
            return self.rad_dms(phi), self.rad_dms(lam), h
             
        elif output == "calc":
             return phi, lam, h 
         
         
    def philamh_XYZ(self, phi, lam, h, inp = 'dec'):
        '''
        funkcja transformujaca wspolrzedne geodezyjne (szerokosc, dlugosc i wysokosc elipsoidalna)
        na wspolrzedne w ukladzie ortokartezjanskim (X,Y,Z)
        
        ---
        parametry:
            phi - szerokosc geodezyjna
            lam - dlugosc geodezyjna
            ^ typ w zaleznosci od inp, mozliwosci:
                    inp = 'dec' -> [float] - stopnie dziesietne (wartosc domyslna)
                    inp = 'dms' -> [tuple] - krotka: (stopnie, minuty, sekundy)
            h [float] - wysokosc 
            
        ---
        returns: 
            X, Y, Z - wspolrzedne ortokartezjanskie 
            
        '''
        if inp == 'dec':
            phi2 = self.degrees_rad(phi)
            lam2 = self.degrees_rad(lam)
            
        elif inp == 'dms':
            phi2 = self.dms_rad(phi)
            lam2 = self.dms_rad(lam)
            
        N = self.obliczenie_N(phi2)
        X1 = (N + h) * m.cos(phi2) * m.cos(lam2)
        Y1 = (N + h) * m.cos(phi2) * m.sin(lam2)
        Z1 = (N*(1-self.e2) + h) * m.sin(phi2)
        
        return X1, Y1, Z1
    
    
    def XYZ_neu(self, X, Y, Z, X0, Y0, Z0):
        """
        

        Parameters
        ----------
        X,Y,Z - [float] współrzędne geocentryczne satelity (te podlegają transformacji)
        X0, Y0, Z0 - [float] współrzędne geocentryczne anteny

        Returns
        -------
       N, E, U - [float] współrzędne satelity w układzie topocentrycznym

        """
       
        phi, lam, _ = [m.radians(coord) for coord in self.XYZ_philamh(X0, Y0, Z0)]
        
        R = np.array([[-m.sin(lam), -m.sin(phi)*m.cos(lam), m.cos(phi)*m.cos(lam)],
                     [m.cos(lam), -m.sin(phi)*m.sin(lam), m.cos(phi)*m.sin(lam)],
                     [0, m.cos(phi), m.sin(phi)]])
        
        xyz_t = np.array([[X-X0],
                         [Y-Y0],
                         [Z-Z0]])
        

        [[E], [N], [U]] = R.T @ xyz_t
        
        return N, E, U
    
    
    def philamh_2000(self, phi, lam, inp = 'dec'):
        '''
        funkcja obliczajaca wspolrzedne w ukladzie pl-2000 na podstawie
        wspolrzednych geodezyjnych punktu (szerokosci i dlugosci elipsoidalnej)
        
        ---
        parametry:
            phi, lam - wspolrzedne geodezyjne punktu, w zaleznosci od
            wartosci inp, moga byc:
                'dec' -> [float] - stopnie dziesietne (domyslne)
                'dms' -> [tuple] - krotka: (stopnie, minuty, sekundy)
                'calc' -> [float] - radiany
        
        ---
        returns: 
            x_2000, y_2000 [float] - wspolrzedne w ukladzie pl-2000 
            
        '''
        if inp == 'dec':
            phi = self.degrees_rad(phi)
            lam = self.degrees_rad(lam)
        elif inp == 'dms':
            phi = self.dms_rad(phi)
            lam = self.dms_rad(lam)
        elif inp == 'calc':
            phi = phi
            lam = lam
        
        if lam > self.degrees_rad(13.5) and lam < self.degrees_rad(16.5):
            lam_0 = self.degrees_rad(15)
            nr_strefy = 5
        elif lam > self.degrees_rad(16.5) and lam < self.degrees_rad(19.5):
            lam_0 = self.degrees_rad(18)
            nr_strefy = 6
        elif lam > self.degrees_rad(19.5) and lam < self.degrees_rad(22.5):
            lam_0 = self.degrees_rad(21)
            nr_strefy = 7
        elif lam > self.degrees_rad(22.5) and lam < self.degrees_rad(25.5):
            lam_0 = self.degrees_rad(24)
            nr_strefy = 8
        
        m2000 = 0.999923
        b2 = (self.a**2)*(1-self.e2)
        e12 = (self.a**2 - b2)/b2
        d_lambda = lam - lam_0
        t = m.tan(phi)
        eta2 = e12 * ((m.cos(phi))**2)
         
        N = self.a/((1 - self.e2 * m.sin(phi)**2)**(1/2))
        sigma = self.dlugosc_luku_poludnika(phi)
            
        x_gk = sigma + (d_lambda**2)/2 * N * m.sin(phi) * m.cos(phi) * (1 + ((d_lambda**2)/12)*((m.cos(phi))**2) * (5 - t**2 + 9 * eta2 + 4 * (eta2**2)) + ((d_lambda**4)/360)*(m.cos(phi))**4 * (61 - 58*(t**2) + t**4 + 270*eta2  - 330 * eta2 * (t**2)))
        y_gk = d_lambda * N * m.cos(phi) * ( 1 + ((d_lambda**2)/6) * m.cos(phi)**2 * (1 - t**2 + eta2)+ ((d_lambda**4)/120)*(m.cos(phi))**4 * (5 - 18*(t**2) + t**4 + 14*eta2  - 58 * eta2 * (t**2)))
        
        x_2000 = x_gk * m2000
        y_2000 = y_gk * m2000 + 500000 + nr_strefy * 1000000
        return x_2000, y_2000
    
    def philamh_1992(self, phi, lam, inp = 'dec'):
        '''
        funkcja obliczajaca wspolrzedne w ukladzie pl-1992 na podstawie
        wspolrzednych geodezyjnych punktu (szerokosci i dlugosci elipsoidalnej)
        
        ---
        parametry:
            phi, lam - wspolrzedne geodezyjne punktu, w zaleznosci od
            wartosci inp, moga byc:
                'dec' -> [float] - stopnie dziesietne (domyslne)
                'dms' -> [tuple] - krotka: (stopnie, minuty, sekundy)
                'calc' -> [float] - radiany
        
        ---
        returns: 
            x_1992, y_1992 [float]- wspolrzedne w ukladzie pl-1992
            
        '''
        if inp == 'dec':
            phi = self.degrees_rad(phi)
            lam = self.degrees_rad(lam)
        elif inp == 'dms':
            phi = self.dms_rad(phi)
            lam = self.dms_rad(lam)
        elif inp == 'calc':
            phi = phi
            lam = lam
        
        lam_0 = self.degrees_rad(19)
        m1992 = 0.9993
        
        b2 = (self.a**2)*(1-self.e2)
        e12 = (self.a**2 - b2)/b2
        d_lambda = lam - lam_0
        t = m.tan(phi)
        eta2 = e12 * ((m.cos(phi))**2)
         
        N = self.a/((1 - self.e2 * m.sin(phi)**2)**(1/2))
        sigma = self.dlugosc_luku_poludnika(phi)
            
        x_gk = sigma + (d_lambda**2)/2 * N * m.sin(phi) * m.cos(phi) * (1 + ((d_lambda**2)/12)*((m.cos(phi))**2) * (5 - t**2 + 9 * eta2 + 4 * (eta2**2)) + ((d_lambda**4)/360)*(m.cos(phi))**4 * (61 - 58*(t**2) + t**4 + 270*eta2  - 330 * eta2 * (t**2)))
        y_gk = d_lambda * N * m.cos(phi) * ( 1 + ((d_lambda**2)/6) * m.cos(phi)**2 * (1 - t**2 + eta2)+ ((d_lambda**4)/120)*(m.cos(phi))**4 * (5 - 18*(t**2) + t**4 + 14*eta2  - 58 * eta2 * (t**2)))
        
        x_1992 = x_gk * m1992 - 5300000
        y_1992 = y_gk * m1992 + 500000
        return x_1992, y_1992
    
        
    def transf_kras(self, xk, yk, zk):
        d11 = 1 - (0.84078048 * 10**(-6))
        d21 = 4.08960007 * 10**(-6)
        d31 = 0.25613864 * 10**(-6)
        d12 = -4.08959962 * 10**(-6)
        d22 = 1 - (0.84078196 * 10**(-6))
        d32 = -1.73888494 * 10**(-6)
        d13 = -0.25614575 * 10**(-6)
        d23 = 1.73888389 * 10**(-6)
        d33 = 1 - (0.84077363 * 10**(-6))
        
        Tx = -33.4297 
        Ty = 146.5746
        Tz = 76.2865
        
        xg = d11 * (xk - Tx) + d12 * (yk - Ty) + d13 * (zk - Tz)
        yg = d21 * (xk - Tx) + d22 * (yk - Ty) + d23 * (zk - Tz)
        zg = d31 * (xk - Tx) + d32 * (yk - Ty) + d33 * (zk - Tz)
        
        return xg, yg, zg
        
    
    
if __name__ == "__main__":
    
    print(sys.argv)

    if sys.argv[2] == 'wgs84':
        model = 'wgs84'
    elif sys.argv[2] == 'grs80':
        model = 'grs80'
    elif sys.argv[2] == 'elipsoida_krasowskiego':
        model = 'elipsoida_krasowskiego'
    elif sys.argv[2] != 'wgs84' and sys.argv[1] != 'grs80' and sys.argv[1] != 'elipsoida_krasowskiego':
        print(f"\n nie rozpoznano modelu {sys.argv[1]}")
        print('\n wprowadz model elipsoidy, ktory jest obslugiwany:\n  > wgs84\n  > grs80\n  > elipsoida_krasowskiego')
        sys.exit(0)
        
    
    header_lines = int(sys.argv[1])
    
    geo = Transformacje(model)
    print(f'\n korzystasz z modelu {model}\n')
    
    input_file_path = sys.argv[-1]
    
    lista_flag = ['--xyz2plh', '--plh2xyz', "--xyz2neu", "--plh2pl2000", "--plh2pl1992", "--XYZkras2XYZgrs80", "--XYZgrs802pl2000", "--XYZgrs802pl1992"]
    
    if len(set(lista_flag) & set(sys.argv)) != 1:
            print(f"\n uzyj JEDNEJ z dostepnych flag:\n > dostepne flagi: {lista_flag}" )
    
    
    elif "--xyz2plh" in sys.argv:
    
        try:
            coords_plh = []
            with open(input_file_path, 'r') as f:
                lines = f.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    x_str, y_str, z_str = line.split(',')
                    x, y, z = float(x_str), float(y_str), float(z_str)
                    p, l, h = geo.XYZ_philamh(x, y, z)
                    coords_plh.append([p, l, h])
                
            with open('result_XYZ_philamh.txt', 'w') as f1:
                f1.write("    phi [deg]      lam [deg]         h [m]\n")
                f1.write(47*'-'+'\n')
                for coords in coords_plh:
                    line = ','.join([f"{coord:^15.6f}" if index!= 2 else f"{coord:^15.3f}" for index, coord in enumerate(coords)])
                    f1.write(line + '\n')
                    
            print("transformacja przebiegla poprawnie\n plik wynikowy dostepny pod nazwa 'result_XYZ_philamh.txt'")
        
        except FileNotFoundError:
            print('sprawdz czy nazwa pliku wejsciowego jest poprawna\n')
                
            
    elif "--plh2xyz" in sys.argv:
        #przy przeliczaniu plh na cos trzeba chyba napisac ze format phi lam do wstawienia to stopnie dziesietne ?
        
        try:
            coords_xyz = []
            with open(input_file_path, 'r') as f2:
                lines = f2.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    phi_str, lam_str, h_str = line.split(',')
                    phi, lam, h = float(phi_str), float(lam_str), float(h_str)
                    x, y, z = geo.philamh_XYZ(phi, lam, h)
                  
                    coords_xyz.append([x, y, z])
                    
            with open('result_philamh_XYZ.txt', 'w') as f4:
                f4.write('     X [m]           Y [m]           Z [m]\n')
                f4.write(47*'-'+'\n')
                for coords in coords_xyz:
                    line = ','.join([f"{coord:^15.3f}" for coord in coords])
                    f4.write(line + '\n')
                    
            print("transformacja przebiegla poprawnie\n plik wynikowy dostepny pod nazwa 'result_philamh_XYZ.txt'")
                    
        except FileNotFoundError:
            print('sprawdz czy nazwa pliku wejsciowego jest poprawna\n')
                
            
    elif "--xyz2neu" in sys.argv:
        
        try:
            coords_neu = []
            with open(input_file_path, 'r') as f5:
                lines = f5.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    x_str, y_str, z_str = line.split(',')
                    x, y, z = float(x_str), float(y_str), float(z_str)
                    X0, Y0, Z0 = [float(coord)for coord in sys.argv[-4:-1]]
                    n, e, u = geo.XYZ_neu(x, y, z, X0, Y0, Z0)
                    coords_neu.append([n, e, u])
                    
            with open('result_XYZ_NEU.txt', 'w') as f6:
                f6.write('     N [m]          E [m]           U [m]\n')
                f6.write(44*'-'+'\n')
                for coords in coords_neu:
                    line = ', '.join([f"{coord:^14.3f}" for coord in coords])
                    f6.write(line + '\n')
                    
            print("transformacja przebiegla poprawnie\n plik wynikowy dostepny pod nazwa 'result_XYZ_NEU.txt'")
                    
        except FileNotFoundError:
            print('sprawdz czy nazwa pliku wejsciowego jest poprawna\n')
    
                
    elif "--plh2pl2000" in sys.argv:
    
        try:
            coords_pl2000 = []
            with open(input_file_path, 'r') as f7:
                lines = f7.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    phi_str, lam_str, h_str = line.split(',')
                    phi, lam, h = float(phi_str), float(lam_str), float(h_str)
                    x_2000, y_2000 = geo.philamh_2000(phi, lam)
                    
                    coords_pl2000.append([x_2000, y_2000])
                    
            with open('result_philamh_pl2000.txt', 'w') as f8:
                f8.write('     X [m]             Y [m]    \n')
                f8.write(33*'-'+'\n')
                for coords in coords_pl2000:
                    line = ' , '.join([f"{coord:^15.3f}" for coord in coords])
                    f8.write(line + '\n')
                    
            print("transformacja przebiegla poprawnie\n plik wynikowy dostepny pod nazwa 'result_philamh_pl2000.txt'")
            
        except FileNotFoundError:
            print('sprawdz czy nazwa pliku wejsciowego jest poprawna\n')
    
            
    elif "--plh2pl1992" in sys.argv:
        
        try:
            coords_pl1992 = []
            with open(input_file_path, 'r') as f9:
                lines = f9.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    phi_str, lam_str, h_str = line.split(',')
                    phi, lam, h = float(phi_str), float(lam_str), float(h_str)
                    x_1992, y_1992 = geo.philamh_1992(phi, lam)
                
                    coords_pl1992.append([x_1992, y_1992])
                    
            with open('result_philamh_pl1992.txt', 'w') as f10:
                f10.write('     X [m]             Y [m]    \n')
                f10.write(32*'-'+'\n')
                for coords in coords_pl1992:
                    line = ' , '.join([f"{coord:^15.3f}"  for coord in coords])
                    f10.write(line + '\n')
            
            print("transformacja przebiegla poprawnie\n plik wynikowy dostepny pod nazwa 'result_philamh_pl1992.txt'")
                    
        except FileNotFoundError:
            print('sprawdz czy nazwa pliku wejsciowego jest poprawna\n')
            
            
    elif "--XYZkras2XYZgrs80" in sys.argv:
        
        try:
            coords_xyzgrs80= []
            with open(input_file_path, 'r') as f11:
                lines = f11.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    x_str, y_str, z_str = line.split(',')
                    x, y, z = float(x_str), float(y_str), float(z_str)
                    xg, yg, zg = geo.transf_kras(x, y, z)
                    coords_xyzgrs80.append([xg, yg, zg])
                    
            with open('result_XYZkras_XYZgrs80.txt', 'w') as f12:
                # plik pomocniczy wiec nie formatujemy, zeby byla wieksza dokladnosc, wiec ma headerlines 0
                for coords in coords_xyzgrs80:
                    line = ','.join([str(coord) for coord in coords])
                    f12.write(line + '\n')
            
            print("transformacja przebiegla poprawnie\n plik wynikowy dostepny pod nazwa 'result_XYZkras_XYZgrs80.txt'")
                    
        except FileNotFoundError:
            print('sprawdz czy nazwa pliku wejsciowego jest poprawna\n')
            
    
    elif "--XYZgrs802pl2000" in sys.argv:
        
        try:
            coords_pl2000_2 = []
            with open(input_file_path, 'r') as f13:
                lines = f13.readlines()
                lines = lines[header_lines:]
                for line in lines:
                    line = line.strip()
                    x_str, y_str, z_str = line.split(',')
                    x, y, z = (float(x_str), float(y_str), float(z_str))
                    phi, lam, h = geo.XYZ_philamh(x, y, z, 'calc')
                    x_2000, y_2000 = geo.philamh_2000(phi, lam, 'calc')
                    
                    coords_pl2000_2.append([x_2000, y_2000])
                    
            with open('result_XYZgrs80_pl2000.txt', 'w') as f14:
                f14.write('     X [m]             Y [m]    \n')
                f14.write(33*'-'+'\n')
                for coords in coords_pl2000_2:
                    line = ' , '.join([f"{coord:^15.3f}" for coord in coords])
                    f14.write(line + '\n')
                    
            print("transformacja przebiegla poprawnie\n plik wynikowy dostepny pod nazwa 'result_XYZgrs80_pl2000.txt'")
            
        except FileNotFoundError:
            print('sprawdz czy nazwa pliku wejsciowego jest poprawna\n')
            
            
    elif "--XYZgrs802pl1992" in sys.argv:
         
         try:
             coords_pl1992_2 = []
             with open(input_file_path, 'r') as f15:
                 lines = f15.readlines()
                 lines = lines[header_lines:] 
                 for line in lines:
                     line = line.strip()
                     x_str, y_str, z_str = line.split(',')
                     x, y, z = (float(x_str), float(y_str), float(z_str))
                     phi, lam, h = geo.XYZ_philamh(x, y, z, 'calc')
                     x_1992, y_1992 = geo.philamh_1992(phi, lam, 'calc')
                     
                     coords_pl1992_2.append([x_1992, y_1992])
                     
             with open('result_XYZgrs80_pl1992.txt', 'w') as f16:
                 f16.write('     X [m]             Y [m]    \n')
                 f16.write(32*'-'+'\n')
                 for coords in coords_pl1992_2:
                     line = ' , '.join([f"{coord:^15.3f}" for coord in coords])
                     f16.write(line + '\n')
                     
             print("transformacja przebiegla poprawnie\n plik wynikowy dostepny pod nazwa 'result_XYZgrs80_pl1992.txt'")
             
         except FileNotFoundError:
             print('sprawdz czy nazwa pliku wejsciowego jest poprawna\n')
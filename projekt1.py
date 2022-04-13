import math
import numpy as np

class transformacje:
    
    def __init__(self, model: str = 'grs80'):
        """
        Parameters
        ----------
        a : [float] - duża półos elipsoidy [m]
        
        b : [float] - mała półos elipsoidy [m]
        
        flattening : [float] - spłaszczenie Ziemii
        
        e2 : [float] - mimosród podniesiony do potęgi 2
        
        -------

        """
        if model == 'wgs84':
            self.a = 6378137.0
            self.b = 6356752.31424528
        elif model == 'grs80':
            self.a = 6378137.0
            self.b = 6356752.31414036 
        else:
            raise NotImplementedError(f"{model} model is not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.e2 = 2 * self.flattening - self.flattening ** 2
          
        
    
    def dms(self, fi):
        '''
        Algorytm zamieniający jednostki kątów ze stopni dziesiętnych
        na stopnie, minuty i sekundy

        Parameters
        ----------
        fi : [float] - dowolny kąt mierzony w radianach 
        
        Returns
        -------
        None.
        
        '''
        d = math.floor(fi)                      #stopnie
        m = math.floor((fi - d) * 60)           #minuty
        s = round((fi - d - (m/60)) * 3600, 5)  #sekundy
        print(d, 'st', m, 'min', s, 'sek')
        
    
            
    def hirvonen(self, X, Y, Z):
        '''
        Algorytm służący do transformacji współrzędnych geocentrycznych ortokartezjańskich (X, Y, Z)
        na współrzędne geodezyjne (fi, lam, h). W wyniku 3-4-krotnej iteracji wyznaczenia wsp. fi można
        przeliczyć współrzędne z dokładnoscią ok 1 cm.
    
        Parameters
        ----------
        X, Y, Z : [float] - współrzędne geocentryczne
        
        Returns
        -------
        fi : [float] - szerokosć geodezyjna [stopnie]
        
        lam : [float] - długosć geodezyjna [stopnie]
        
        h : [float] - wysokosć elipsoidalna [m]
        
        '''
        
        #XYZ2blh
        r = math.sqrt(X**2 + Y**2)
        fi_n = math.atan(Z/(r * (1 - self.e2)))
        eps = 0.000001/3600 * math.pi/180
        fi = fi_n * 2
        
        while np.abs (fi_n - fi) > eps:
            fi = fi_n
            N = self.a / math.sqrt(1 - self.e2 * math.sin(fi_n)**2)
            h = r / math.cos(fi_n) - N
            fi_n = math.atan(Z/(r * (1 - self.e2 * (N/(N + h)))))
        
        lam = math.atan(Y/X)
        N = self.a / math.sqrt(1 - self.e2 * math.sin(fi_n)**2)
        h = r / math.cos(fi_n) - N
        
        fi_n = fi_n * 180/math.pi
        lam = lam * 180/math.pi 
        
        return fi_n, lam, h
    
    
    
    def blh2xyz(self, fi, lam, h):
        '''
        Zamiana współrzędnych geodezyjnych (fi, lam, h) na współrzędne geocentryczne 
        ortokartezjańskie (X, Y, Z). Algorytm odwrtotny do Hirvonena.
        
        Parameters
        ----------
        fi : [float] - szerokosć geodezyjna [stopnie]
        
        lam : [float] - długosć geodezyjna [stopnie]
        
        h : [float] - wysokosć elipsoidalna [m]

        Returns
        -------
        X, Y, Z : [float] - współrzędne geocentryczne

        '''
        
        fi = fi * math.pi/180
        lam = lam * math.pi/180
        
        N = self.a / math.sqrt(1 - self.e2 * (math.sin(fi))**2)
        X = (N + h) * math.cos(fi) * math.cos(lam)
        Y = (N + h) * math.cos(fi) * math.sin(lam)
        Z = (N * (1 - self.e2) + h) * math.sin(fi)
        
        return X, Y, Z
    
    
    
    def xyz2neu(self, X1, Y1, Z1, X2, Y2, Z2):
        """
        Zamiana współrzędnych geocentrycznych ortokartezjańskich (X, Y, Z), na współrzędne
        topocentryczne (n, e, u)

        Parameters
        ----------
        X1, Y1, Z1 : [float] - współrzędne ortokartezjańskie pierwszego punktu
        
        X2, Y2, Z2 : [float] - współrzędne ortokartezjańskie drugiego punktu

        Returns
        -------
        n, e, u : [float] - współrzędne topocentryczne

        """
        
        delta_X = X1 - X2
        delta_Y = Y1 - Y2
        delta_Z = Z1 - Z2
        
        fi, lam, h = self.hirvonen(X1, Y1, Z1)
        fi = fi * math.pi/180
        lam = lam * math.pi/180
        
        R = np.array([[-np.sin(fi) * np.cos(lam), -np.sin(lam), np.cos(fi) * np.cos(lam)],
                      [-np.sin(fi) * np.sin(lam),  np.cos(lam), np.cos(fi) * np.sin(lam)],
                      [ np.cos(fi), 0, np.sin(fi)]])
        
        d = np.matrix([delta_X, delta_Y, delta_Z])
        d = d.T
        neu = R * d
        
        n = float(neu[0])
        e = float(neu[1])
        u = float(neu[2])
        
        return n, e, u
    
    
    
    
    def uklad1992(self, fi, lam):
        """
        Zamiana współrzędnych geodezyjnych (fi, lam, h) na współrzędne płaskie 
        polskiego układu 1992 (x, y)

        Parameters
        ----------
        fi : [float] - szerokoć geograficzna [stopnie]
        
        lam : [float] - długosć geograficzna [stopnie]

        Returns
        -------
        x92, y92 : [float] - współrzędne płaskie w układzie 1992

        """
        
        fi = fi * math.pi/180
        lam = lam * math.pi/180
        
        N = self.a/(math.sqrt(1 - self.e2 * np.sin(fi)**2))
        t = np.tan(fi)
        n2 = self.e2 * np.cos(lam)**2
        lam0 = math.radians(19)
        l = lam - lam0
        m0 = 0.9993
        
        A0 = 1 - (self.e2/4) - (3*(self.e2**2))/64 - (5 * (self.e2**3))/256
        A2 = 3/8 * (self.e2 + ((self.e2**2)/4) + ((15*self.e2**3)/128))
        A4 = 15/256 * (self.e2**2 + (3 * (self.e2**3))/4)
        A6 = (35 * (self.e2**3))/3072
        
        sigma = self.a * ((A0 * fi) - (A2 * np.sin(2 * fi)) + (A4 * np.sin(4 * fi)) - (A6 * np.sin(6 * fi)))
        
        x = sigma + ((l**2)/2) * (N * np.sin(fi) * np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + (9 * n2) + (4 * n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58 * (t**2)) + (t**4) + (270 * n2) - (330 * n2 *(t**2))))
        y = l * (N * np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1 - (t**2) + n2)) + (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14 * n2) - (58 * n2 * (t**2))))

        x92 = round(x * m0 - 5300000, 4)
        y92 = round(y * m0 + 500000, 4)   
        
        return x92, y92
    
    
    def uklad2000(self, fi, lam):
        """
        Zamiana współrzędnych geodezyjnych (fi, lam, h) na współrzędne płaskie 
        polskiego układu 2000 (x, y)

        Parameters
        ----------
        fi : [float] - szerokoć geograficzna [stopnie]
        
        lam : [float] - długosć geograficzna [stopnie]

        Returns
        -------
        x00, y00 : [float] - współrzędne płaskie w układzie 2000

        """
        
        fi = fi * math.pi/180
        lam = lam * math.pi/180
        
        m0 = 0.999923
        N = self.a/(math.sqrt((1 - self.e2) * np.sin(fi)**2))
        t = np.tan(fi)
        n2 = self.e2 * np.cos(lam)**2
        lam = math.degrees(lam)
        
        if lam > 13.5 and lam < 16.5:
            s = 5
            lam0 = 15
        elif lam > 16.5 and lam < 19.5:
            s = 6
            lam0 = 18
        elif lam > 19.5 and lam < 22.5:
            s = 7
            lam0 = 21
        elif lam > 22.5 and lam < 25.5:
            s = 8
            lam0 = 24
        elif lam < 13.5 and lam > 25.5:
            print('punktu nie da się przeniesć do układu 2000')
            
        lam = math.radians(lam)
        lam0 = math.radians(lam0)
        l = lam - lam0
        
        A0 = 1 - (self.e2/4) - (3 * (self.e2**2))/64 - (5 * (self.e2**3))/256
        A2 = 3/8 * (self.e2 + ((self.e2**2)/4) + ((15 * (self.e2**3))/128))
        A4 = 15/256 * (self.e2**2 + (3 * (self.e2**3))/4)
        A6 = (35 * (self.e2**3))/3072
        
        
        sigma = self.a * ((A0 * fi) - (A2 * np.sin(2 * fi)) + (A4 * np.sin(4 * fi)) - (A6 * np.sin(6 * fi)))
        
        x = sigma + ((l**2)/2) * (N * np.sin(fi) * np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + (9 * n2) + (4 * (n2**2))) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58 * (t**2)) + (t**4) + (270 * n2) - (330 * n2 * (t**2))))
        y = l * (N * np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1 - (t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14 * n2) - (58 * n2 * (t**2))))

        x00 = round(x * m0, 4)
        y00 = round(y * m0 + (s * 1000000) + 500000, 4)   
        
        return x00, y00 
    
    
    def katy_odl(self, X1, Y1, Z1, X2, Y2, Z2):
        """
        Algorytm wyznaczający kąt azymutu, kąt elewacji, oraz odległosci 2D i 3D
        dla układu topocentrycznego (n, e, u)

        Parameters
        ----------
        X1, Y1, Z1 : [float] - współrzędne ortokartezjańskie pierwszego punktu
        
        X2, Y2, Z2 : [float] - współrzędne ortokartezjańskie drugiego punktu

        Returns
        -------
        Az : [float] - azymut [stopnie]
        
        alfa : [float] - kąt elewacji [stopnie]
        
        dist2D : [float] - odległosc 2D [m]
        
        dist3D : [float] - odległosć 3D [m]

        """
        
        n, e, u = self.xyz2neu(X1, Y1, Z1, X2, Y2, Z2)
        
        
        Az = math.atan2(e, n)
        
        if Az < 0:
            Az = Az + 2 * math.pi
        else:
            pass
        
        dist2D = math.sqrt(e**2 + u**2)
        dist3D = math.sqrt(n**2 + e**2 + u**2)
        alfa = math.atan2(u, dist2D)
        
        Az = Az * 180/math.pi
        alfa = alfa * 180/math.pi
        
        return Az, alfa, dist2D, dist3D
    
    
        
    
if __name__ == "__main__":
    
    # utworzenie obiektu
    test = transformacje(model = "grs80")
    
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    
    # transformacje
    fi, lam, h = test.hirvonen(X, Y, Z)
    print('geodezyjne:', 'fi =', fi, 'lam =', lam, 'h =', h)
    
    X, Y, Z = test.blh2xyz(fi, lam, h)
    print('geocentryczne:', 'X =', round(X, 5), 'Y =', round(Y, 5), 'Z =', round(Z, 5))
    
    x92, y92 = test.uklad1992(fi, lam)
    print("współrzędne plaskie układu 1992:", 'x92 =', round(x92, 5), 'y92 =', round(y92, 5)) 
    
    x00, y00 = test.uklad2000(fi, lam)
    print("współrzędne plaskie układu 2000:", 'x00 =', round(x00, 5), 'y00 =', round(y00, 5))
    
    n, e, u = test.xyz2neu(X, Y, Z, X + 456, Y + 10, Z + 8962)
    print('współrzędne topocentryczne:', 'n =', round(n, 5), 'e =', round(e, 5), 'u =', round(u, 5))
    
    Az, alfa, dist2D, dist3D = test.katy_odl(X, Y, Z, X + 456, Y + 10, Z + 8962)
    print('Az =', Az)
    print('alfa =', alfa)
    print('odległosć 2D =', round(dist2D, 5))
    print('odległosć 3D =', round(dist3D, 5))
    
    print('fi =')
    fi = test.dms(fi)
    
    print('lam =')
    lam = test.dms(lam)
    
    print('Azymut =')
    Az = test.dms(Az)
    
    print('kąt elewacji =')
    alfa = test.dms(alfa)
    
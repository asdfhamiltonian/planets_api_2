'''
Planetary positions calculator
Programmed with algorithms/equations from "How to Compute Planetary Positions"
by Paul Schlyter http://www.stjarnhimlen.se/comp/ppcomp.html#20
'''

from math import copysign, floor, pi, atan2, sqrt, cos, sin, acos, asin
from functools import partial
import datetime
from collections import OrderedDict


def degrees(rads):
    """converts radians to degrees"""
    return rads / pi * 180


def radians(degrees):
    """converts degrees to radians"""
    return degrees / 180 * pi


def fnday(y, m, d, h):
    """calculates epoch day from calendar date"""
    a = 367 * y - 7 * (y + (m + 9) // 12) // 4 + 275 * m // 9
    a += d - 730530 + h / 24
    return a


def fnipart(x):
    """returns the whole number portion of a real number, positive or negative"""
    return copysign(1, x) * floor(abs(x))


def fnrange(angle):
    """converts angle in radians to the range 0 to 2pi"""
    newangle = angle % (2 * pi)
    return newangle


def fnatan2(y, x):
    """atan2 function"""
    a = atan2(y, x)
    if a < 0:
        a = a + 2 * pi
    return a


def timeadjust(utcdis, utc):
    """returns local time, given UTC displacement and UTC time"""
    x = utcdis + utc
    while x < 0:
        x = x + 24
    if x > 24:
        x = x - 24
    return x


def decimaltominsecs(x):
    """decimal of an hour time to minutes and seconds"""
    x = x % 1
    mins = floor(x * 60)
    x = (x * 60) % 1
    secs = floor(x * 60)
    return (mins, secs)


class Planets(object):
    localLat = 0
    localLong = 0

    def __init__(self, llat=42, llon=-122):
        self.localLat = llat
        self.localLong = llon

    def calcecl(self, day):
        ecl = radians(23.4393 - 3.563E-07 * day)
        return ecl

    def calcposition(self, N, i, w, a, e, m, Epoch=2000.0, day=0):
        bigE1 = 1
        bigE0 = 0
        bvar = 1

        bigE = m + e * sin(m) * (1.0 + e * cos(m))

        if e > 0.05:
            while abs(bvar) > 1E-5:
                ratio = (bigE - e * sin(bigE) - m) / (1 - e * cos(bigE))
                bigE1 = bigE - ratio
                bvar = degrees(bigE1) - degrees(bigE)
                bigE = bigE1
        xv = a * (cos(bigE) - e)
        yv = a * (sqrt(1 - e * e) * sin(bigE))

        v = fnatan2(yv, xv)
        r = sqrt(xv * xv + yv * yv)

        xh = r * (cos(N) * cos(v + w) - sin(N) * sin(v + w) * cos(i))
        yh = r * (sin(N) * cos(v + w) + cos(N) * sin(v + w) * cos(i))
        zh = r * (sin(v + w) * sin(i))
        lonecl = fnatan2(yh, xh)
        latecl = fnatan2(zh, sqrt(xh * xh + yh * yh))
        lonCorr = 3.82394E-5 * (365.2422 * (Epoch - 2000.0) - day)

        return [xh, yh, zh, lonecl, latecl, v, r, Epoch, lonCorr]

    def calcsun(self, day):
        sunArray = self.Sun(day)
        Nsun = sunArray[0]
        isun = sunArray[1]
        wsun = sunArray[2]
        asun = sunArray[3]
        esun = sunArray[4]
        Msun = sunArray[5]
        ecl = radians(23.4393 - 3.563E-07 * day)
        positionArray = self.calcposition(Nsun, isun, wsun, asun, esun, Msun)
        vsun = positionArray[5]
        rs = positionArray[6]
        lonsun = fnrange(vsun + wsun)
        meanlonsun = fnrange(Msun + wsun)
        '''meanlonsun is different from lonsun'''

        xs = rs * cos(lonsun)
        ys = rs * sin(lonsun)
        zs = 0

        """E = Msun + esun * sin(Msun) * (1.0 + esun * cos(Msun))
        xv = cos(E) - esun
        yv = sqrt(1.0 - esun*esun) * sin(E)

        vsun = fnatan2(yv, xv)
        rs = sqrt(xv*xv + yv*yv)

        #lonsun = Msun + wsun
        lonsun = fnrange(vsun + wsun)

        xs = rs * cos(lonsun)
        ys = rs * sin(lonsun)"""

        '''sun rectangular coordinates'''
        xe = xs
        ye = ys * cos(ecl)
        ze = ys * sin(ecl)

        RA = fnatan2(ye, xe)
        Dec = fnatan2(ze, sqrt(xe * xe + ye * ye))

        return [rs, lonsun, vsun, xs, ys, xe, ye, ze, RA, Dec, meanlonsun]

    def calcGeocentric(self, lonecl, latecl, ecl, r, xs, ys):
        xh = r * cos(lonecl) * cos(latecl)
        yh = r * sin(lonecl) * cos(latecl)
        zh = r * sin(latecl)

        xg = xh + xs
        yg = yh + ys
        zg = zh

        xe = xg
        ye = yg * cos(ecl) - zg * sin(ecl)
        ze = yg * sin(ecl) + zg * cos(ecl)

        ra = fnatan2(ye, xe)
        dec = fnatan2(ze, sqrt(xe * xe + ye * ye))

        rg = sqrt(xe * xe + ye * ye + ze * ze)

        return [xe, ye, ze, ra, dec, rg]

    def calcTopocentric(self, h, mins, second, day, meanlonsun,
                        r, ra, dec, lat, lon):
        h = h + mins / 60 + second / 3600
        lat1 = radians(lat)
        gmsto = 12 * (meanlonsun + pi) / pi
        gmst = gmsto + h
        lst = gmst + self.localLong / 15
        lst = lst % 24
        ha = lst - degrees(ra) / 15
        while (ha > 24) or (ha < -24):
            if ha > 12:
                ha = ha - 24
            elif ha < -12:
                ha = ha + 24
        ha = radians(ha * 15)
        x = cos(ha) * cos(dec)
        y = sin(ha) * cos(dec)
        z = sin(dec)

        xhor = x * sin(lat1) - z * cos(lat1)
        yhor = y
        zhor = x * cos(lat1) + z * sin(lat1)

        az = atan2(yhor, xhor) + pi
        alt = asin(zhor)

        ppar = radians(8.794 / 3600) / r
        alttopoc = alt - ppar * cos(alt)

        return [az, alttopoc, alt, ppar]

    def calcGeoTopo(self, h, mins, second, d, lat, lon,
                    lonplanet, latplanet, r):
        '''
        Takes in parameters for time, position, ecliptic position. Returns
        azimuth, altitude, right ascension, declination.
        '''
        sunArray = self.calcsun(d)  # sun's position
        '''Returns the following array:
        [rs, lonsun, vsun, xs, ys, xe, ye, ze, RA, Dec, meanlonsun]'''
        lonsun = sunArray[1]
        meanlonsun = sunArray[10]
        xs = sunArray[3]
        ys = sunArray[4]
        ecl = self.calcecl(d)

        positionArray = self.calcGeocentric(lonplanet, latplanet,
                                            ecl, r, xs, ys)
        '''returns [xe, ye, ze, ra, dec, rg]'''
        xe = positionArray[0]
        ye = positionArray[1]
        ze = positionArray[2]
        ra = positionArray[3]
        dec = positionArray[4]
        re = positionArray[5]

        topArray = self.calcTopocentric(h, mins, second, d, meanlonsun,
                                        re, ra, dec, lat, lon)
        return [degrees(topArray[0]), degrees(topArray[1]),
                degrees(ra) / 15, dec]

    def calcElements(self, N, i, w, a, e, M):
        w1 = N + w  # longitude of perihelion
        L = M + w1  # mean longitude
        q = a * (1 - e)  # perihelion distance
        Q = a * (1 + e)  # aphelion distance
        P = a ^ 1.5  # orbital period (years if in AU)
        '''T = Epoch_of_M - (M(deg)/360)) / P = time of perihelion'''

    def timeNow(self):
        y = datetime.datetime.utcnow().year
        m = datetime.datetime.utcnow().month
        d = datetime.datetime.utcnow().day
        h = datetime.datetime.utcnow().hour
        mins = datetime.datetime.utcnow().minute
        second = datetime.datetime.utcnow().second
        return y, m, d, h, mins, second

    def dateNow(self):
        y = datetime.datetime.utcnow().year
        m = datetime.datetime.utcnow().month
        d = datetime.datetime.utcnow().day
        return y, m, d

    def Sun(self, day):
        Nsun = 0.0
        isun = 0.0
        wsun = fnrange(radians(282.9404 + 4.70935E-5 * day))
        asun = 1.000000  # (AU)
        esun = 0.016709 - 1.151E-9 * day
        Msun = fnrange(radians(356.0470 + 0.9856002585 * day))
        return [Nsun, isun, wsun, asun, esun, Msun]

    def Moon(self, day):
        Nm = fnrange(radians(125.1228 - 0.0529538083 * day))
        im = radians(5.1454)
        wm = fnrange(radians(318.0634 + 0.1643573223 * day))
        am = 60.266  # earth radii
        em = 0.054900
        Mm = fnrange(radians(115.3654 + 13.0649929509 * day))
        return [Nm, im, wm, am, em, Mm]

    def Mercury(self, day):
        Nmer = fnrange(radians(48.3313 + 3.24587E-5 * day))
        imer = fnrange(radians(7.0047 + 5.00E-8 * day))
        wmer = fnrange(radians(29.1241 + 1.01444E-5 * day))
        amer = 0.387098  # (AU)
        emer = 0.205635 + 5.59E-10 * day
        Mmer = fnrange(radians(168.6562 + 4.0923344368 * day))
        return [Nmer, imer, wmer, amer, emer, Mmer]

    def Mars(self, day):
        Nmar = fnrange(radians(49.5574 + 2.11081E-5 * day))
        imar = fnrange(radians(1.8497 - 1.78E-8 * day))
        wmar = fnrange(radians(286.5016 + 2.92961E-5 * day))
        amar = 1.523688  # (AU)
        emar = 0.093405 + 2.516E-9 * day
        Mmar = fnrange(radians(18.6021 + 0.5240207766 * day))
        return [Nmar, imar, wmar, amar, emar, Mmar]

    def Venus(self, day):
        Nven = fnrange(radians(76.6799 + 2.46590E-5 * day))
        iven = fnrange(radians(3.3946 + 2.75E-8 * day))
        wven = fnrange(radians(54.8910 + 1.38374E-5 * day))
        aven = 0.723330  # (AU)
        even = 0.006773 - 1.302E-9 * day
        Mven = fnrange(radians(48.0052 + 1.6021302244 * day))
        return [Nven, iven, wven, aven, even, Mven]

    def Jupiter(self, day):
        Nj = fnrange(radians(100.4542 + 2.76854E-5 * day))
        ij = fnrange(radians(1.3030 - 1.557E-7 * day))
        wj = fnrange(radians(273.8777 + 1.64505E-5 * day))
        aj = 5.20256  # (AU)
        ej = 0.048498 + 4.469E-9 * day
        Mj = fnrange(radians(19.8950 + 0.0830853001 * day))
        return [Nj, ij, wj, aj, ej, Mj]

    def Saturn(self, day):
        Nsat = fnrange(radians(113.6634 + 2.38980E-5 * day))
        isat = fnrange(radians(2.4886 - 1.081E-7 * day))
        wsat = fnrange(radians(339.3939 + 2.97661E-5 * day))
        asat = 9.55475  # (AU)
        esat = 0.055546 - 9.499E-9 * day
        Msat = fnrange(radians(316.9670 + 0.0334442282 * day))
        return [Nsat, isat, wsat, asat, esat, Msat]

    def Uranus(self, day):
        Nu = fnrange(radians(74.0005 + 1.3978E-5 * day))
        iu = fnrange(radians(0.7733 - 1.9E-8 * day))
        wu = fnrange(radians(96.6612 + 3.0565E-5 * day))
        au = 19.18171 - 1.55E-8 * day  # (AU)
        eu = 0.047318 + 7.45E-9 * day
        Mu = fnrange(radians(142.5905 + 0.011725806 * day))
        return [Nu, iu, wu, au, eu, Mu]

    def Neptune(self, day):
        Nnep = fnrange(radians(131.7806 + 3.0173E-5 * day))
        inep = fnrange(radians(1.7700 - 2.55E-7 * day))
        wnep = fnrange(radians(272.8461 - 6.027E-6 * day))
        anep = 30.05826 + 3.313E-8 * day  # (AU)
        enep = 0.008606 + 2.15E-9 * day
        Mnep = fnrange(radians(260.2471 + 0.005995147 * day))
        return [Nnep, inep, wnep, anep, enep, Mnep]

    def calcMars(self, day):
        mArray = self.Mars(day)

        Nm = mArray[0]
        im = mArray[1]
        wm = mArray[2]
        am = mArray[3]
        em = mArray[4]
        Mm = mArray[5]

        calcArray = self.calcposition(Nm, im, wm, am, em, Mm)
        lon = calcArray[3]
        lat = calcArray[4]
        r = calcArray[6]
        return [lon, lat, r]

    def calcNeptune(self, day):
        nepArray = self.Neptune(day)

        N = nepArray[0]
        i = nepArray[1]
        w = nepArray[2]
        a = nepArray[3]
        e = nepArray[4]
        M = nepArray[5]

        calcArray = self.calcposition(N, i, w, a, e, M)
        lon = calcArray[3]
        lat = calcArray[4]
        r = calcArray[6]
        return [lon, lat, r]

    def calcMercury(self, day):
        merArray = self.Mercury(day)
        N = merArray[0]
        i = merArray[1]
        w = merArray[2]
        a = merArray[3]
        e = merArray[4]
        M = merArray[5]

        calcArray = self.calcposition(N, i, w, a, e, M)
        lon = calcArray[3]
        lat = calcArray[4]
        r = calcArray[6]
        return [lon, lat, r]

    def calcVenus(self, day):
        vArray = self.Venus(day)

        Nv = vArray[0]
        iv = vArray[1]
        wv = vArray[2]
        av = vArray[3]
        ev = vArray[4]
        Mv = vArray[5]

        calcArray = self.calcposition(Nv, iv, wv, av, ev, Mv)
        lon = calcArray[3]
        lat = calcArray[4]
        r = calcArray[6]
        return [lon, lat, r]

    def calcJupiter(self, day):
        jArray = self.Jupiter(day)
        satArray = self.Saturn(day)
        uArray = self.Uranus(day)
        Mj = jArray[5]
        Ms = satArray[5]
        Mu = uArray[5]

        Nj = jArray[0]
        ij = jArray[1]
        wj = jArray[2]
        aj = jArray[3]
        ej = jArray[4]
        '''Mj = jArray[5]'''

        calcArray = self.calcposition(Nj, ij, wj, aj, ej, Mj)
        lon = calcArray[3]
        lat = calcArray[4]
        r = calcArray[6]
        dlon = 0

        '''these are in degrees, added at end'''

        dlon -= 0.332 * sin(2 * Mj - 5 * Ms - radians(67.6))
        dlon -= 0.056 * sin(2 * Mj - 2 * Ms + radians(21))
        dlon += 0.042 * sin(3 * Mj - 5 * Ms + radians(21))
        dlon -= 0.036 * sin(Mj - 2 * Ms)
        dlon += 0.022 * cos(Mj - Ms)
        dlon += 0.023 * sin(2 * Mj - 3 * Ms + radians(52))
        dlon -= 0.016 * sin(Mj - 5 * Ms + radians(69))

        lon = radians(dlon) + lon
        return [lon, lat, r]

    def calcSaturn(self, day):
        jArray = self.Jupiter(day)
        satArray = self.Saturn(day)
        uArray = self.Uranus(day)
        Mj = jArray[5]
        Ms = satArray[5]
        Mu = uArray[5]

        Nsat = satArray[0]
        isat = satArray[1]
        wsat = satArray[2]
        asat = satArray[3]
        esat = satArray[4]
        '''Mj = jArray[5]'''

        calcArray = self.calcposition(Nsat, isat, wsat, asat, esat, Ms)
        lon = calcArray[3]
        lat = calcArray[4]
        r = calcArray[6]
        dlon = 0
        dlat = 0

        dlon += 0.812 * sin(2 * Mj - 5 * Ms - radians(67.6))
        dlon -= 0.229 * cos(2 * Mj - 4 * Ms - radians(2))
        dlon += 0.119 * sin(Mj - 2 * Ms - radians(3))
        dlon += 0.046 * sin(2 * Mj - 6 * Ms - radians(69))
        dlon += 0.014 * sin(Mj - 3 * Ms + radians(32))

        dlat -= 0.020 * cos(2 * Mj - 4 * Ms - radians(2))
        dlat += 0.018 * sin(2 * Mj - 6 * Ms - radians(49))

        lon = radians(dlon) + lon
        lat = radians(dlat) + lat

        return [lon, lat, r]

    def calcUranus(self, day):
        jArray = self.Jupiter(day)
        satArray = self.Saturn(day)
        uArray = self.Uranus(day)
        Mj = jArray[5]
        Ms = satArray[5]
        Mu = uArray[5]

        Nu = uArray[0]
        iu = uArray[1]
        wu = uArray[2]
        au = uArray[3]
        eu = uArray[4]
        '''Mj = jArray[5]'''

        calcArray = self.calcposition(Nu, iu, wu, au, eu, Mu)
        lon = calcArray[3]
        lat = calcArray[4]
        r = calcArray[6]
        dlon = 0

        dlon += 0.040 * sin(Ms - 2 * Mu + radians(6))
        dlon += 0.035 * sin(Ms - 3 * Mu + radians(33))
        dlon -= 0.015 * sin(Mj - Mu + radians(20))

        lon = radians(dlon) + lon

        return [lon, lat, r]

    def calcMercuryNow(self):
        y, m, d, h, mins, second = self.timeNow()
        return self.calcMercuryDate(y, m, d, h, mins, second)

    def calcVenusNow(self):
        y, m, d, h, mins, second = self.timeNow()
        return self.calcVenusDate(y, m, d, h, mins, second)

    def calcJupiterNow(self):
        y, m, d, h, mins, second = self.timeNow()
        return self.calcJupiterDate(y, m, d, h, mins, second)

    def calcSaturnNow(self):
        y, m, d, h, mins, second = self.timeNow()
        return self.calcSaturnDate(y, m, d, h, mins, second)

    def calcUranusNow(self):
        y, m, d, h, mins, second = self.timeNow()
        return self.calcUranusDate(y, m, d, h, mins, second)

    def calcNeptuneNow(self):
        y, m, d, h, mins, second = self.timeNow()
        return self.calcNeptuneDate(y, m, d, h, mins, second)

    def calcMarsNow(self):
        y, m, d, h, mins, second = self.timeNow()
        return self.calcMarsDate(y, m, d, h, mins, second)

    def calcMercuryDate(self, y, m, d, h, mins, second):
        lat = self.localLat
        lon = self.localLong
        d = fnday(y, m, d, h)
        d += mins / (60 * 24) + second / (3600 * 24)
        merArray = self.calcMercury(d)
        ecl = self.calcecl(d)
        lonmer = merArray[0]
        latmer = merArray[1]
        rmer = merArray[2]
        return self.calcGeoTopo(h, mins, second, d, lat, lon,
                                lonmer, latmer, rmer)

    def calcVenusDate(self, y, m, d, h, mins, second):
        lat = self.localLat
        lon = self.localLong
        d = fnday(y, m, d, h)
        d += mins / (60 * 24) + second / (3600 * 24)
        vArray = self.calcVenus(d)
        ecl = self.calcecl(d)
        lonv = vArray[0]
        latv = vArray[1]
        rv = vArray[2]

        return self.calcGeoTopo(h, mins, second, d, lat, lon,
                                lonv, latv, rv)

    def calcMarsDate(self, y, m, d, h, mins, second):
        lat = self.localLat
        lon = self.localLong
        d = fnday(y, m, d, h)
        d += mins / (60 * 24) + second / (3600 * 24)
        mArray = self.calcMars(d)
        ecl = self.calcecl(d)
        lonm = mArray[0]
        latm = mArray[1]
        rm = mArray[2]

        return self.calcGeoTopo(h, mins, second, d, lat, lon,
                                lonm, latm, rm)

    def calcJupiterDate(self, y, m, d, h, mins, second):
        lat = self.localLat
        lon = self.localLong
        d = fnday(y, m, d, h)
        d += mins / (60 * 24) + second / (3600 * 24)
        jArray = self.calcJupiter(d)
        ecl = self.calcecl(d)
        lonj = jArray[0]
        latj = jArray[1]
        rj = jArray[2]

        return self.calcGeoTopo(h, mins, second, d, lat, lon,
                                lonj, latj, rj)

    def calcSaturnDate(self, y, m, d, h, mins, second):
        lat = self.localLat
        lon = self.localLong
        d = fnday(y, m, d, h)
        d += mins / (60 * 24) + second / (3600 * 24)
        satArray = self.calcSaturn(d)
        ecl = self.calcecl(d)
        lonsat = satArray[0]
        latsat = satArray[1]
        rsat = satArray[2]

        return self.calcGeoTopo(h, mins, second, d, lat, lon,
                                lonsat, latsat, rsat)

    def calcUranusDate(self, y, m, d, h, mins, second):
        lat = self.localLat
        lon = self.localLong
        d = fnday(y, m, d, h)
        d += mins / (60 * 24) + second / (3600 * 24)
        uArray = self.calcUranus(d)
        ecl = self.calcecl(d)
        lonu = uArray[0]
        latu = uArray[1]
        ru = uArray[2]

        return self.calcGeoTopo(h, mins, second, d, lat, lon,
                                lonu, latu, ru)

    def calcNeptuneDate(self, y, m, d, h, mins, second):
        lat = self.localLat
        lon = self.localLong
        d = fnday(y, m, d, h)
        d += mins / (60 * 24) + second / (3600 * 24)
        nArray = self.calcNeptune(d)
        ecl = self.calcecl(d)
        lonn = nArray[0]
        latn = nArray[1]
        rn = nArray[2]

        return self.calcGeoTopo(h, mins, second, d, lat, lon,
                                lonn, latn, rn)

    def calcPlutoThatIsNotAPlanetNow(self):
        lat = self.localLat
        lon = self.localLong
        y = datetime.datetime.utcnow().year
        m = datetime.datetime.utcnow().month
        d = datetime.datetime.utcnow().day
        h = datetime.datetime.utcnow().hour
        mins = datetime.datetime.utcnow().minute
        second = datetime.datetime.utcnow().second
        d = fnday(y, m, d, h)
        d += mins / (60 * 24) + second / (3600 * 24)
        '''
        now to do some calculations for pluto's position based on numerical
        integration based calculations for pluto's orbit, following
        Paul Schlyter's example
        '''
        S = radians(50.03 + 0.033459652 * d)
        P = radians(238.95 + 0.003968789 * d)

        lonecl = radians(238.9508 + 0.00400703 * d)
        lonadd = -19.799 * sin(P) + 19.848 * cos(P)
        lonadd += 0.897 * sin(2 * P) - 4.956 * cos(2 * P)
        lonadd += 0.610 * sin(3 * P) + 1.211 * cos(3 * P)
        lonadd += -0.341 * sin(4 * P) - 0.190 * cos(4 * P)
        lonadd += 0.128 * sin(5 * P) - 0.034 * cos(5 * P)
        lonadd += -0.038 * sin(6 * P) + 0.031 * cos(6 * P)
        lonadd += 0.020 * sin(S - P) - 0.010 * cos(S - P)
        lonecl += radians(lonadd)

        latecl = radians(-3.9082)
        latadd = -5.453 * sin(P) - 14.975 * cos(P)
        latadd += 3.527 * sin(2 * P) + 1.673 * cos(2 * P)
        latadd += -1.051 * sin(3 * P) + 0.328 * cos(3 * P)
        latadd += 0.179 * sin(4 * P) - 0.292 * cos(4 * P)
        latadd += 0.019 * sin(5 * P) + 0.100 * cos(5 * P)
        latadd += -0.031 * sin(6 * P) - 0.026 * cos(6 * P)
        latadd += 0.011 * cos(S - P)
        latecl += radians(latadd)

        r = 40.72
        radd = 6.68 * sin(P) + 6.90 * cos(P)
        radd += -1.18 * sin(2 * P) - 0.03 * cos(2 * P)
        radd += +0.15 * sin(3 * P) - 0.14 * cos(3 * P)
        r += radians(radd)

        return self.calcGeoTopo(h, mins, second, d, lat, lon,
                                lonecl, latecl, r)

    def calcPlutoDate(self, y, m, d, h, mins, second):
        lat = self.localLat
        lon = self.localLong
        d = fnday(y, m, d, h)
        d += mins / (60 * 24) + second / (3600 * 24)
        '''
        now to do some calculations for pluto's position from numerical
        integration derived parameters for pluto's orbit, following
        Paul Schlyter's example
        '''
        S = radians(50.03 + 0.033459652 * d)
        P = radians(238.95 + 0.003968789 * d)

        lonecl = radians(238.9508 + 0.00400703 * d)
        lonadd = -19.799 * sin(P) + 19.848 * cos(P)
        lonadd += 0.897 * sin(2 * P) - 4.956 * cos(2 * P)
        lonadd += 0.610 * sin(3 * P) + 1.211 * cos(3 * P)
        lonadd += -0.341 * sin(4 * P) - 0.190 * cos(4 * P)
        lonadd += 0.128 * sin(5 * P) - 0.034 * cos(5 * P)
        lonadd += -0.038 * sin(6 * P) + 0.031 * cos(6 * P)
        lonadd += 0.020 * sin(S - P) - 0.010 * cos(S - P)
        lonecl += radians(lonadd)

        latecl = radians(-3.9082)
        latadd = -5.453 * sin(P) - 14.975 * cos(P)
        latadd += 3.527 * sin(2 * P) + 1.673 * cos(2 * P)
        latadd += -1.051 * sin(3 * P) + 0.328 * cos(3 * P)
        latadd += 0.179 * sin(4 * P) - 0.292 * cos(4 * P)
        latadd += 0.019 * sin(5 * P) + 0.100 * cos(5 * P)
        latadd += -0.031 * sin(6 * P) - 0.026 * cos(6 * P)
        latadd += 0.011 * cos(S - P)
        latecl += radians(latadd)

        r = 40.72
        radd = 6.68 * sin(P) + 6.90 * cos(P)
        radd += -1.18 * sin(2 * P) - 0.03 * cos(2 * P)
        radd += +0.15 * sin(3 * P) - 0.14 * cos(3 * P)
        r += radians(radd)

        return self.calcGeoTopo(h, mins, second, d, lat, lon,
                                lonecl, latecl, r)

    def localTime(self, t, utcdis):
        t += utcdis
        t = t % 24
        return t

    def localTimeArray(self, array, utcdis):
        local = partial(self.localTime, utcdis=utcdis)
        return list(map(local, array))

    def calcRiseSet(self, calcPlanetDate, utcdis, y, m, d):
        lat = self.localLat
        lon = self.localLong
        h = 12 - utcdis  # calculate for noon local time
        epochday = fnday(y, m, d, h)  # not same as day of month
        sunArray = self.calcsun(epochday)
        planetArray = calcPlanetDate(y, m, d, h, 0, 0)
        ra = planetArray[2]
        decl = planetArray[3]
        Ls = sunArray[10]  # mean longitude of the sun
        gmst0 = 12 * (Ls + pi) / (pi)  # gmst0
        b = radians(-0.583)
        '''degrees below horizon for set/rise, accounting for refraction'''
        UTPlanetInSouth = ((ra * 15) - (gmst0 * 15) - lon) / 15.0
        if UTPlanetInSouth < 0:
            UTPlanetInSouth += 24
        cosLHA = (sin(b) - sin(radians(lat)) * sin(decl))
        cosLHA /= (cos(radians(lat)) * cos(decl))
        if cosLHA < -1.0:
            raise ValueError('Planet always above altitude limit')
        elif cosLHA > 1.0:
            raise ValueError('Planet always below altitude limit')
        else:
            LHA = degrees(acos(cosLHA)) / 15.04107
            planetrise = UTPlanetInSouth - LHA
            planetset = UTPlanetInSouth + LHA
            utcArray = [planetrise, planetset,
                        UTPlanetInSouth, LHA]  # UTC time
            return self.localTimeArray(utcArray, utcdis)  # local time

    def calcMercuryRiseSet(self, utcdis):
        y, m, d = self.dateNow()
        return self.calcRiseSet(self.calcMercuryDate, utcdis, y, m, d)

    def calcVenusRiseSet(self, utcdis):
        y, m, d = self.dateNow()
        return self.calcRiseSet(self.calcVenusDate, utcdis, y, m, d)

    def calcMarsRiseSet(self, utcdis):
        y, m, d = self.dateNow()
        return self.calcRiseSet(self.calcMarsDate, utcdis, y, m, d)

    def calcJupiterRiseSet(self, utcdis):
        y, m, d = self.dateNow()
        return self.calcRiseSet(self.calcJupiterDate, utcdis, y, m, d)

    def calcSaturnRiseSet(self, utcdis):
        y, m, d = self.dateNow()
        return self.calcRiseSet(self.calcSaturnDate, utcdis, y, m, d)

    def calcUranusRiseSet(self, utcdis):
        y, m, d = self.dateNow()
        return self.calcRiseSet(self.calcUranusDate, utcdis, y, m, d)

    def calcNeptuneRiseSet(self, utcdis):
        y, m, d = self.dateNow()
        return self.calcRiseSet(self.calcNeptuneDate, utcdis, y, m, d)

    def calcPlutoRiseSet(self, utcdis):
        y, m, d = self.dateNow()
        return self.calcRiseSet(self.calcPlutoDate, utcdis, y, m, d)

    def calcMercuryRiseSetDate(self, utcdis, y, m, d):
        return self.calcRiseSet(self.calcMercuryDate, utcdis, y, m, d)

    def calcVenusRiseSetDate(self, utcdis, y, m, d):
        return self.calcRiseSet(self.calcVenusDate, utcdis, y, m, d)

    def calcMarsRiseSetDate(self, utcdis, y, m, d):
        return self.calcRiseSet(self.calcMarsDate, utcdis, y, m, d)

    def calcJupiterRiseSetDate(self, utcdis, y, m, d):
        return self.calcRiseSet(self.calcJupiterDate, utcdis, y, m, d)

    def calcSaturnRiseSetDate(self, utcdis, y, m, d):
        return self.calcRiseSet(self.calcSaturnDate, utcdis, y, m, d)

    def calcUranusRiseSetDate(self, utcdis, y, m, d):
        return self.calcRiseSet(self.calcUranusDate, utcdis, y, m, d)

    def calcNeptuneRiseSetDate(self, utcdis, y, m, d):
        return self.calcRiseSet(self.calcNeptuneDate, utcdis, y, m, d)

    def calcPlutoRiseSetDate(self, utcdis, y, m, d):
        return self.calcRiseSet(self.calcPlutoDate, utcdis, y, m, d)

    def allRiseSet(self, utcdis):
        risesets = [self.calcMercuryRiseSet(utcdis),
                    self.calcVenusRiseSet(utcdis),
                    self.calcMarsRiseSet(utcdis),
                    self.calcJupiterRiseSet(utcdis),
                    self.calcSaturnRiseSet(utcdis),
                    self.calcUranusRiseSet(utcdis),
                    self.calcNeptuneRiseSet(utcdis),
                    self.calcPlutoRiseSet(utcdis)]
        planetnames = ["mercury", "venus", "mars", "jupiter",
                       "saturn", "uranus", "neptune", "pluto"]
        riseSetString = ""
        for num in range(0, 8):
            data = risesets[num]
            a = planetnames[num]
            a += ": rise {}, set {}. \n".format(data[0], data[1])
            riseSetString += a

        return riseSetString

    def dictRiseSet(self, utcdis):
        risesets = [self.calcMercuryRiseSet(utcdis),
                    self.calcVenusRiseSet(utcdis),
                    self.calcMarsRiseSet(utcdis),
                    self.calcJupiterRiseSet(utcdis),
                    self.calcSaturnRiseSet(utcdis),
                    self.calcUranusRiseSet(utcdis),
                    self.calcNeptuneRiseSet(utcdis),
                    self.calcPlutoRiseSet(utcdis)]
        planetnames = ["mercury", "venus", "mars", "jupiter",
                       "saturn", "uranus", "neptune", "pluto"]
        planetDict = OrderedDict()
        for num in range(0, 8):
            data = risesets[num]
            subdict = {}
            subdict["rise"] = data[0]
            subdict["set"] = data[1]
            planetDict[planetnames[num]] = subdict
        return planetDict

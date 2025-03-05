import datetime
import json
import logging
import numpy as np

logging.basicConfig(level=logging.INFO)

class mpcorb_file():
    """
    Read and write MPCORB files ussing the format stated in https://www.minorplanetcenter.net/iau/info/MPOrbitFormat.html on march 4, 2025 reproduced below:

    Export Format for Minor-Planet Orbits
    This document describes the format used for both unperturbed and perturbed orbits of minor planets, as used in the Extended Computer Service and in the Minor Planet Ephemeris Service.

    Orbital elements for minor planets are heliocentric.

    The column headed `F77' indicates the Fortran 77/90/95/2003/2008 format specifier that should be used to read the specified value.

    Columns   F77    Use

        1 -   7  a7     Number or provisional designation
                        (in packed form)
        9 -  13  f5.2   Absolute magnitude, H
    15 -  19  f5.2   Slope parameter, G

    21 -  25  a5     Epoch (in packed form, .0 TT)
    27 -  35  f9.5   Mean anomaly at the epoch, in degrees

    38 -  46  f9.5   Argument of perihelion, J2000.0 (degrees)
    49 -  57  f9.5   Longitude of the ascending node, J2000.0
                        (degrees)
    60 -  68  f9.5   Inclination to the ecliptic, J2000.0 (degrees)

    71 -  79  f9.7   Orbital eccentricity
    81 -  91  f11.8  Mean daily motion (degrees per day)
    93 - 103  f11.7  Semimajor axis (AU)

    106        i1     Uncertainty parameter, U
            or a1     If this column contains `E' it indicates
                        that the orbital eccentricity was assumed.
                        For one-opposition orbits this column can
                        also contain `D' if a double (or multiple)
                        designation is involved or `F' if an e-assumed
                        double (or multiple) designation is involved.

    108 - 116  a9     Reference
    118 - 122  i5     Number of observations
    124 - 126  i3     Number of oppositions

        For multiple-opposition orbits:
        128 - 131  i4     Year of first observation
        132        a1     '-'
        133 - 136  i4     Year of last observation

        For single-opposition orbits:
        128 - 131  i4     Arc length (days)
        133 - 136  a4     'days'

    138 - 141  f4.2   r.m.s residual (")
    143 - 145  a3     Coarse indicator of perturbers
                        (blank if unperturbed one-opposition object)
    147 - 149  a3     Precise indicator of perturbers
                        (blank if unperturbed one-opposition object)
    151 - 160  a10    Computer name

    There may sometimonth be additional information beyond column 160
    as follows:

    162 - 165  z4.4   4-hexdigit flags

                        The bottom 6 bits (bits 0 to 5) are used to encode
                        a value representing the orbit type (other
                        values are undefined):

                        Value
                            1  Atira
                            2  Aten
                            3  Apollo
                            4  Amor
                            5  Object with q < 1.665 AU
                            6  Hungaria
                            7  Unused or internal MPC use only
                            8  Hilda
                            9  Jupiter Trojan
                        10  Distant object

                        Additional information is conveyed by
                        adding in the following bit values:

                Bit  Value
                    6     64  Unused or internal MPC use only
                    7    128  Unused or internal MPC use only
                    8    256  Unused or internal MPC use only
                    9    512  Unused or internal MPC use only
                    10   1024  Unused or internal MPC use only
                    11   2048  Object is NEO
                    12   4096  Object is 1-km (or larger) NEO
                    13   8192  1-opposition object seen at
                            earlier opposition
                    14  16384  Critical list numbered object
                    15  32768  Object is PHA

                        Note that the orbit classification is
                        based on cuts in osculating element
                        space and is not 100% reliable.

                        Note also that certain of the flags
                        are for internal MPC use and are
                        not documented.

    167 - 194  a      Readable designation

    195 - 202  i8     Date of last observation included in
                        orbit solution (YYYYMMDD format)

    """

    def __init__(self,file=None):
        self.bodies=None
        if not file is None:
            self.read(file)     
        
    # Función para convertir el formato comprimido de la época a fecha juliana
    def compressed_epoch_to_datetime(self,epoch:str)->datetime.datetime:
        """
        Convert compressed epoch to python datetime following the below rules:

        Dates of the form YYYYMMDD may be packed into five characters to conserve space.

        The first two digits of the year are packed into a single character in column 1 (I = 18, J = 19, K = 20). 
        Columns 2-3 contain the last two digits of the year.
        Column 4 contains the month and column 5 contains the day, coded as detailed below:

        Month     Day      Character         Day      Character
                            in Col 4 or 5              in Col 4 or 5
        Jan.       1           1             17           H
        Feb.       2           2             18           I
        Mar.       3           3             19           J
        Apr.       4           4             20           K
        May        5           5             21           L
        June       6           6             22           M
        July       7           7             23           N
        Aug.       8           8             24           O
        Sept.      9           9             25           P
        Oct.      10           A             26           Q
        Nov.      11           B             27           R
        Dec.      12           C             28           S
                  13           D             29           T
                  14           E             30           U
                  15           F             31           V
                  16           G

        Examples:

        1996 Jan. 1    = J9611
        1996 Jan. 10   = J961A
        1996 Sept.30   = J969U
        1996 Oct. 1    = J96A1
        2001 Oct. 22   = K01AM

        This system can be extended to dates with non-integral days. The decimal fraction of the day is simply appended to the five characters defined above.

        Examples:

        1998 Jan. 18.73     = J981I73
        2001 Oct. 22.138303 = K01AM138303
        """
        year_letter = {'I':'18','J':'19','K':'20'}
        month_map = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9,
                'A': 10, 'B': 11, 'C': 12}
        day_map={
                '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9,
                'A': 10, 'B': 11, 'C': 12,'D': 13, 'E': 14, 'F': 15,
                'G': 16, 'H': 17, 'I': 18,'J': 19, 'K': 20, 'L': 21,
                'M': 22, 'N': 23, 'O': 24,'P': 25, 'Q': 26, 'R': 27,
                'S': 28, 'T': 29, 'U': 30,'V': 31
                }

        year_letter_epoch = epoch[0]
        month_epoch = epoch[3]
        day_epoch = epoch[4]
        
        if year_letter.get(year_letter_epoch, 0)==0:
            raise ValueError(f"Invalid epoch format: {epoch}")
        year = f'{year_letter.get(year_letter_epoch, 0)}{epoch[1:3]}'
        month = month_map.get(month_epoch, 0)
        day = day_map.get(day_epoch, 0)
        if month == 0 or day == 0:
            raise ValueError(f"Invalid epoch format: {epoch}")
        date_str=f"{year}-{month:02d}-{day:02d}"
        date=datetime.datetime.strptime(date_str,'%Y-%m-%d')
        return date

    def body_to_numeric(self,body:dict)->dict:
        '''
        converted to numeric values instead of its string representatio
        '''
        changed_values=dict()
        for k,v in body.items():
            if k in ['G','H','M','Peri','Node','i','e','n','a']:
                try:
                    changed_values[k]=float(v)
                except ValueError:
                    changed_values[k]=np.nan
        #changed_values['epoch_datetime']=self.compressed_epoch_to_datetime(body['epoch'])
        body.update(changed_values)
        return body
    
    def to_numeric(self)->list:
        '''
        Return a list of dictinaries with variables converted to numeric values instead of its string representation
        '''
        result=list()
        for body in self.bodies:
                result.append(self.body_to_numeric(body))
        return result
    
    def parse_line(self,l:str)->dict:
        """
        Parse one line an return a dict with all the variables fullfilled.
        """
        # Orbital elements
        line=' '+l   #padding to sync index with mpcorb description
        nombre = line[1:8].strip()
        G=line[9:14].strip()
        H=line[15:20].strip()
        epoch = line[21:26].strip()  # Época en formato comprimido
        M = line[27:36].strip()  # Anomalía meday (grados)
        omega = line[38:47].strip()  # Argumento del perihelio (grados)
        Omega = line[49:58].strip()  # Longitud del nodo ascendente (grados)
        i = line[60:69].strip()  # Inclinación (grados)
        e = line[71:80].strip()  # Excentricidad
        n = line[81:92].strip()  # Movimiento medio (grados/día)
        a = line[93:104].strip()  # Semieje mayor (AU)    
        rest=line[105:]
        
        body={
            'nombre': nombre,
            'G':G,
            'H':H,
            'a': a,
            'e': e,
            'i': i,
            'n':n,
            'Peri': omega,
            'Node': Omega,
            'M': M,
            'epoch': epoch,
            'rest':rest
            }     
        return body

    def make_line(self,body:dict)->str:
        """
        Compose one line with the body data 
        """        
        # Ceres data used to dim line
        ceres='00001    3.34  0.15 K2555 188.70269   73.27343   80.25221   10.58780  0.0794013  0.21424651   2.7660512  0 E2024-V47  7330 125 1801-2024 0.80 M-v 30k MPCLINUX   4000      (1) Ceres              20241101'
        line=[' ' for x in range(len(ceres))]
        line[1:8]=body['nombre'].ljust(7)
        line[9:14]=body['G'].rjust(5)
        line[15:20]=body['H'].rjust(5)
        line[21:26]=body['epoch'].rjust(5)
        line[27:36]=body['M'].rjust(9)
        line[38:47]=body['Peri'].rjust(9)
        line[49:58]=body['Node'].rjust(9)
        line[60:69]=body['i'].rjust(9)
        line[71:80]=body['e'].rjust(9)
        line[81:92]=body['n'].rjust(11)
        line[93:104]=body['a'].rjust(11)
        line[105:]=body['rest'].replace('\n','')
        line=line[1:]
        result=''.join(line)
        return result     

    def read(self,filename:str)->list:
        '''
        Read the MPCORB.DAT file.
        '''
        with open(filename,'r') as fd:
            lines=fd.readlines()
        #skip header if any (all text above '---')    
        start_line=[ i for i,l in enumerate(lines) if '---' in l]

        if len(start_line)>0:
            lines=lines[start_line[-1]+1:]

        #load all bodies
        bodies=list()
        for l in lines:
            if l.startswith('#') or len(l.strip())<1:  # Ignore empty lines or comments
                continue
            body=self.parse_line(l)
            bodies.append(body)
        self.bodies=bodies  #save classwise to caching when called by other fn
        return bodies
    
    def write(self,filename:str):
         '''
         Write a file formated as MPCORB with the bodies data
         '''
         if self.bodies==None:
            return False
         with open(filename,'w') as fd:
            for body in self.bodies:
                fd.write(self.make_line(body))
                fd.write('\n')
            return True
         
    def write_json(self,filename):
        with open(filename, 'w') as f:
            json.dump(self.bodies, f,indent=2)

import datetime
import json
import logging
import numpy as np
from tqdm import tqdm
from typing import Generator

logging.basicConfig(level=logging.INFO)

def json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""
    if isinstance(obj, (datetime.datetime, datetime.date)):
        return obj.isoformat()
    raise TypeError("Type %s not serializable" % type(obj))


def hook(obj):
    """ tqdm hook for json.loads """
    #TBD
    return obj

# convenience fn
def add_asteroids_to_rebound(simulation, bodies=None):
    """
    Add asteroids to a REBOUND simulation.

    Example::

        import rebound
        import mpcorbfile
        import numpy as np

        sim=rebound.Simulation()
        rebound.data.add_solar_system(sim)

        mpcorb = 'MPCORB_TEST.DAT'
        f = mpcorbfile.mpcorb_file(mpcorb)
        mpcorbfile.add_asteroids_to_rebound(sim)

    """

    for body in bodies:
        # append to simulation
        simulation.add(
            m=0,  # Masa del cuerpo (0 para asteroides)
            a=body["a"],
            e=body["e"],
            inc=np.radians(body["i"]),
            omega=np.radians(body["Peri"]),
            Omega=np.radians(body["Node"]),
            M=np.radians(body["M"]),
            date=body["epochJD"],
            hash=body["Name"],
        )


def set_elliptical_body_elements(eliptical_body, body):
    """
    Set orbital elements of eliptical_body for futher calculation
    using pyephem::

        pyephem parameters:

        _inc        — Inclination (°)
        _Om         — Longitude of ascending node (°)
        _om         — Argument of perihelion (°)
        _a          — Mean distance from sun (AU)
        _M          — Mean anomaly from the perihelion (°)
        _epoch_M    — Date for measurement _M
        _size       — Angular size (arcseconds at 1 AU)
        _e          — Eccentricity
        _epoch      — Epoch for _inc, _Om, and _om
        _H, _G      — Parameters for the H/G magnitude model
        _g, _k      — Parameters for the g/k magnitude model

        eliptical_body._H = body["H"]
        eliptical_body._G = body["G"]
        eliptical_body._a = body["a"]
        eliptical_body._M = body["M"]
        eliptical_body._om = body["Peri"]
        eliptical_body._Om = body["Node"]
        eliptical_body._inc = body["i"]
        eliptical_body._e = body["e"]
        eliptical_body._epoch = body["Epoch"].strftime("%Y/%m/%d %H:%M:%S")
        eliptical_body._epoch_M = body["Epoch"].strftime("%Y/%m/%d %H:%M:%S")

        
    """
    eliptical_body._H = body["H"]
    eliptical_body._G = body["G"]
    eliptical_body._a = body["a"]
    eliptical_body._M = body["M"]
    eliptical_body._om = body["Peri"]
    eliptical_body._Om = body["Node"]
    eliptical_body._inc = body["i"]
    eliptical_body._e = body["e"]
    eliptical_body._epoch = body["Epoch"].strftime("%Y/%m/%d %H:%M:%S")
    eliptical_body._epoch_M = body["Epoch"].strftime("%Y/%m/%d %H:%M:%S")
    return eliptical_body


class mpcorb_file:
    """
    Read and write MPCORB files ussing the format stated in 
    https://www.minorplanetcenter.net/iau/info/MPOrbitFormat.html on march 4, 2025
    
    format::

        Export Format for Minor-Planet Orbits
        This document describes the format used for both unperturbed and perturbed orbits of minor planets,
        as used in the Extended Computer Service and in the Minor Planet Ephemeris Service.

        Orbital elements for minor planets are heliocentric.

        The column headed `F77' indicates the Fortran 77/90/95/2003/2008 format specifier
        that should be used to read the specified value.

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

    def __init__(self, file=None):
        self.bodies = list()
        # This fields dont throw exceptions if missing
        self.optional_fields = [
            "G",
            "H",
            "U",
            "Ref",
            "Num_obs",
            "Num_opps",
            "Arc_length",
            "rms",
            "Perturbers",
            "Perturbers_2",
            "Computer",
            "Hex_flags",
            "Last_obs",
            "Number",
            "Name",
        ]
        # if type is not present field is consider string
        self.format_dict = {
            "packed_designation": {"from": 1, "to": 8, "ljust": True},
            "H": {"from": 9, "to": 14, "type": float, "ljust": False, "format": "5.2f"},
            "G": {
                "from": 15,
                "to": 20,
                "type": float,
                "ljust": False,
                "format": "5.2f",
            },
            "Epoch": {
                "from": 21,
                "to": 26,
                "type": self.compressed_epoch,
                "ljust": False,
                "format": "",
            },
            "M": {
                "from": 27,
                "to": 36,
                "type": float,
                "ljust": False,
                "format": "9.5f",
            },
            "Peri": {
                "from": 38,
                "to": 47,
                "type": float,
                "ljust": False,
                "format": "9.5f",
            },
            "Node": {
                "from": 49,
                "to": 58,
                "type": float,
                "ljust": False,
                "format": "9.5f",
            },
            "i": {
                "from": 60,
                "to": 69,
                "type": float,
                "ljust": False,
                "format": "9.5f",
            },
            "e": {
                "from": 71,
                "to": 80,
                "type": float,
                "ljust": False,
                "format": "9.7f",
            },
            "n": {
                "from": 81,
                "to": 92,
                "type": float,
                "ljust": False,
                "format": "11.8f",
            },
            "a": {
                "from": 93,
                "to": 104,
                "type": float,
                "ljust": False,
                "format": "11.7f",
            },
            "U": {"from": 106, "to": 107, "ljust": False, "format": ""},
            "Ref": {"from": 108, "to": 117, "ljust": False, "format": ""},
            "Num_obs": {
                "from": 118,
                "to": 123,
                "type": int,
                "ljust": False,
                "format": "5d",
            },
            "Num_opps": {
                "from": 124,
                "to": 127,
                "type": int,
                "ljust": False,
                "format": "3d",
            },
            "Arc_length": {"from": 128, "to": 137, "ljust": False, "format": ""},
            "rms": {
                "from": 138,
                "to": 142,
                "type": float,
                "ljust": True,
                "format": "4.2f",
            },
            "Perturbers": {"from": 143, "to": 146, "ljust": False, "format": ""},
            "Perturbers_2": {"from": 147, "to": 150, "ljust": True, "format": ""},
            "Computer": {"from": 151, "to": 161, "ljust": True, "format": ""},
            "Hex_flags": {
                "from": 162,
                "to": 166,
                "type": "hex",
                "ljust": False,
                "format": "04X",
            },
            "Number": {"from": 167, "to": 175, "ljust": False, "format": ""},
            "Name": {"from": 176, "to": 194, "ljust": True, "format": ""},
            "Last_obs": {"from": 195, "to": 203, "ljust": False, "format": ""},
        }
        if file is not None:
            self.read(file)

    def datetime_compressed_epoch(self, epoch: datetime.datetime) -> str:
        year_letter_map = {"18": "I", "19": "J", "20": "K"}
        day_map = "123456789ABCDEFGHIJKLMNOPQRSTUV"
        year_str = str(epoch.year)
        year_letter = year_letter_map[year_str[0:2]]
        result = f"{year_letter}{year_str[2:4]}{epoch.month:1X}{day_map[epoch.day - 1]}"
        return result

    def compressed_epoch(self, epoch: str) -> datetime.datetime:
        return self.compressed_epoch_to_datetime(epoch)

    # Función para convertir el formato comprimido de la época a fecha juliana
    def compressed_epoch_to_datetime(self, epoch: str) -> datetime.datetime:
        """
        Convert compressed epoch to python datetime following the below rules::

            Dates of the form YYYYMMDD may be packed into five characters to conserve space.

            The first two digits of the year are packed into a single character in 
            column 1 (I = 18, J = 19, K = 20).
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

            This system can be extended to dates with non-integral days. The decimal fraction of 
            the day is simply appended to the five characters defined above.

            Examples:

            1998 Jan. 18.73     = J981I73
            2001 Oct. 22.138303 = K01AM138303

        """
        year_letter = {"I": "18", "J": "19", "K": "20"}
        month_map = {
            "1": 1,
            "2": 2,
            "3": 3,
            "4": 4,
            "5": 5,
            "6": 6,
            "7": 7,
            "8": 8,
            "9": 9,
            "A": 10,
            "B": 11,
            "C": 12,
        }
        day_map = {
            "1": 1,
            "2": 2,
            "3": 3,
            "4": 4,
            "5": 5,
            "6": 6,
            "7": 7,
            "8": 8,
            "9": 9,
            "A": 10,
            "B": 11,
            "C": 12,
            "D": 13,
            "E": 14,
            "F": 15,
            "G": 16,
            "H": 17,
            "I": 18,
            "J": 19,
            "K": 20,
            "L": 21,
            "M": 22,
            "N": 23,
            "O": 24,
            "P": 25,
            "Q": 26,
            "R": 27,
            "S": 28,
            "T": 29,
            "U": 30,
            "V": 31,
        }

        year_letter_epoch = epoch[0]
        month_epoch = epoch[3]
        day_epoch = epoch[4]

        if year_letter.get(year_letter_epoch, 0) == 0:
            raise ValueError(f"Invalid epoch format: {epoch}")
        year = f"{year_letter.get(year_letter_epoch, 0)}{epoch[1:3]}"
        month = month_map.get(month_epoch, 0)
        day = day_map.get(day_epoch, 0)
        if month == 0 or day == 0:
            raise ValueError(f"Invalid epoch format: {epoch}")
        date_str = f"{year}-{month:02d}-{day:02d}"
        date = datetime.datetime.strptime(date_str, "%Y-%m-%d")
        return date

    def __add_calculate_fields(self, body: dict) -> dict:
        """
        Add some calculated fields to the body dict
        """
        newbody = body.copy()
        newbody["epochJD"] = self.datetime_to_julian_date(newbody["Epoch"])

        # Add calculate new fields
        # print(newbody['packed_designation'])
        newbody["designation"] = self.expand_packed_designation(
            newbody["packed_designation"]
        )
        # newbody["discover_date"] = self.date_from_packed_designation(newbody['packed_designation'])
        # Better use designation. Packed designation losses his date meaning when asteroid get numbered
        #'name' field has discovery date meaning while it is provisional (no given name).
        newbody["discover_date"] = self.date_from_designation(newbody["Name"])
        newbody["orbit_type"] = self.orbit_type(
            newbody["a"], newbody["e"], newbody["i"]
        )
        return newbody

    def datetime_to_julian_date(self, my_date: datetime.datetime) -> float:
        """
        Convert a datetime to julian date
        """
        return my_date.toordinal() + 1721424.5

    def add(self, body_dict: dict):
        """
        Add new body from a dict.
        """
        # TODO check keys
        _body_dict = self.__add_calculate_fields(body_dict)
        self.bodies.append(_body_dict)

    def __parse_line(self, line: str) -> dict:
        """
        Parse one line an return a dict with all the variables fullfilled.
        """
        # line = " " + l  # padding to sync index with mpcorb description
        body = dict()
        for k, v in self.format_dict.items():
            if "type" in v:
                try:
                    if v["type"] == "hex":
                        body[k] = int(line[v["from"] - 1 : v["to"] - 1], 16)
                    else:
                        body[k] = v["type"](line[v["from"] - 1 : v["to"] - 1])
                except ValueError:
                    body[k] = np.nan
            else:
                body[k] = line[v["from"] - 1 : v["to"] - 1].strip()
        body = self.__add_calculate_fields(body)
        return body

    def __make_line(self, body: dict) -> str:
        """
        Compose one line with the body data
        """

        # Ceres data used to dim line
        ceres = "00001    3.34  0.15 K2555 188.70269   73.27343   80.25221   10.58780  0.0794013  0.21424651   2.7660512  0 E2024-V47  7330 125 1801-2024 0.80 M-v 30k MPCLINUX   4000      (1) Ceres              20241101"
        line = [" " for x in range(len(ceres))]
        if True:
            for k, v in self.format_dict.items():
                if k not in body:
                    if k in self.optional_fields:
                        # Fill with default values
                        if "type" not in v:
                            body[k] = ""
                        elif v["type"] == "hex":
                            body[k] = 0x0
                        elif v["type"] is float:
                            body[k] = v["type"](np.nan)
                        else:
                            body[k] = 0
                    else:
                        raise Exception(f"Required field: {k} not in body:{body}")

                if "format" in v:
                    if "type" in v:
                        if v["type"] == self.compressed_epoch:
                            txt = self.datetime_compressed_epoch(body[k])
                        elif np.isnan(body[k]):
                            txt = ""
                        else:
                            txt = f"{body[k]:{v['format']}}"
                    else:
                        txt = f"{body[k]:{v['format']}}"
                else:
                    txt = body[k]
                if v["ljust"]:
                    text = txt.ljust(v["to"] - v["from"])
                else:
                    text = txt.rjust(v["to"] - v["from"])
                line[v["from"] - 1 : v["to"] - 1] = text
        # except:
        #    print('Fail to typer',k,v,body[k],type(body[k]))
        return "".join(line)

    def read(self, filename: str) -> list:
        """
        Read the MPCORB.DAT file.
        """
        bodies = []
        with open(filename, "r") as fd:
            lines = fd.readlines()
        # skip header if any (all text above '---')
        start_line = [i for i, line in enumerate(lines) if "---" in line]

        if len(start_line) > 0:
            lines = lines[start_line[-1] + 1 :]

        # load all bodies
        bodies = list()
        for line in tqdm(
            lines, colour="green", unit=" bodies", desc="reading", unit_scale=True
        ):
            if (
                line.startswith("#") or len(line.strip()) < 1
            ):  # Ignore empty lines or comments
                continue
            body = self.__parse_line(line)
            bodies.append(body)
        self.bodies = bodies  # save classwise to caching when called by other fn
        self.colnames = list(bodies[0].keys())
        return bodies

    def read_json(self, filename: str) -> list:
        """
        read json files https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz
        """
        with open(filename, "r") as f:
            self.bodies = json.load(f, object_hook=hook)
            # TO BE DONE
            for body in tqdm(self.bodies):
                if "Name" not in body:
                    body["Name"] = body["Principal_desig"]
                if "Number" in body:
                    body["packed_designation"] = self.pack_designation(body["Number"])
                elif "Name" in body:
                    body["packed_designation"] = self.pack_designation(body["Name"])
                else:
                    body["packed_designation"] = ""
                body["Epoch"] = datetime.date.fromordinal(
                    int(body["Epoch"] - 1721424.5)
                )  # From julian date
                if "Arc_years" in body:
                    body["Arc_length"] = body["Arc_years"]
                elif "Arc_length" in body:
                    body["Arc_length"] = f"{body['Arc_length']} days"
                    pass
                else:
                    body["Arc_length"] = ""
                body["Last_obs"] = body["Last_obs"].replace("-", "")
                body["Hex_flags"] = int(body["Hex_flags"], 16)
        return self.bodies

    def write(self, filename: str, header: str = "") -> bool:
        """
        Write a file formated as MPCORB with the bodies data
        """
        if self.bodies is None:
            return False
        Note = "\n                               Create with mpcorbfile python library/utility. See: https://github.com/nachoplus/mpcorbfile\n"
        colnames = "Des'n     H     G   Epoch     M        Peri.      Node       Incl.       e            n           a        Reference #Obs #Opp    Arc    rms  Perts   Computer"
        with open(filename, "w") as fd:
            fd.write(f"{header}\n")
            fd.write(f"{Note}\n")
            fd.write(f"{colnames}\n")
            fd.write("".join(["-" for x in range(len(colnames) + 2)]))
            fd.write("\n")
            for body in tqdm(
                self.bodies,
                colour="blue",
                unit=" bodies",
                desc="writting",
                unit_scale=True,
            ):
                fd.write(self.__make_line(body))
                fd.write("\n")
            return True

    def write_json(self, filename: str):
        """Write json file compatible with https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz" files"""
        with open(filename, "w") as f:
            json.dump(self.bodies, f, indent=2, default=json_serial)

    def json(self):
        return json.dumps(self.bodies, indent=2, default=json_serial)

    def get_chunks(self, n: int) -> list:
        """return a n list of lists with len(list)/n bodies each"""
        N = int(np.ceil(len(self.bodies) / n))
        return self.__group(self.bodies, N)

    # Internal fn
    def __group(self, lst: list, n: int) -> Generator[list, None, None]:
        for i in range(0, len(lst), n):
            val = lst[i : i + n]
            yield val

    def __hex2dec(self, letter: chr) -> str:
        """
        Convert 0..F hex digit to 0..15 decimal
        """
        try:
            int(letter)
            return letter
        except ValueError:
            if letter.isupper():
                return str(ord(letter) - ord("A") + 10)
            if letter.islower():
                return str(ord(letter) - ord("a") + 36)

    def pack_designation(self, designation: str) -> str:
        """
        Create packed designation from designation following format: https://www.minorplanetcenter.net/iau/info/PackedDes.html
        """
        year_letter_map = {"18": "I", "19": "J", "20": "K"}
        base62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
        isdigit = str.isdigit
        if designation[0] == "(" and designation[-1] == ")":
            # Numered designation
            number = int(designation[1:-1])
            if number < 100000:
                result = f"{number:05d}"
            elif number < 620000:
                first_two = int(number / 10000)
                modu = number % 10000
                result = f"{base62[first_two]}{modu:04d}"
            else:
                number = number - 620000
                bit4, rest = divmod(number, 62**3)
                bit3, rest = divmod(rest, 62**2)
                bit2, rest = divmod(rest, 62**1)
                bit1, rest = divmod(rest, 1)
                # print(bit1,bit2,bit3,bit4,rest)
                result = f"~{base62[bit4]}{base62[bit3]}{base62[bit2]}{base62[bit1]}"
            return result
        elif isdigit(designation[0:4]) and designation[4] == " ":
            # Provisional designation
            # Survey
            if designation[5:8] in ["P-L", "T-1", "T-2", "T-3"]:
                return f"{designation[5:8].replace('-', '')}{designation[0:4]}"
            # Regular
            year = f"{year_letter_map[designation[0:2]]}{designation[2:4]}"
            halfmonth = designation[5:6]
            letter = designation[6:7]
            if len(designation) > 7:
                number = int(designation[7:])
                # print(year,halfmonth,letter,number)
                if int(number) > 99:
                    # first two digits
                    first_two = designation[7:9]
                    last = designation[9]
                    number_txt = f"{base62[int(first_two)]}{int(last)}"
                    result = f"{year}{halfmonth}{number_txt}{letter}"
                else:
                    result = f"{year}{halfmonth}{number:02d}{letter}"
            else:
                result = f"{year}{halfmonth}00{letter}"
            return result

    def expand_packed_designation(self, packed: str) -> str:
        """
        Convert the packed designation format to formal designation following format: https://www.minorplanetcenter.net/iau/info/PackedDes.html
        """

        isdigit = str.isdigit
        desig = ""
        try:
            packed = packed.strip()
        except ValueError:
            print("ValueError: Input is not convertable to string.")

        if isdigit(packed):
            desig = packed.lstrip("0")  # ex: 00123
        elif (
            not isdigit(packed[0]) and packed[0] != "~"
        ):  # ex: A7659 = 107659 but not ~0000
            if isdigit(packed[1:]):  # ex: A7659
                desig = self.__hex2dec(packed[0]) + packed[1:]

            elif isdigit(packed[1:3]):  # ex: J98SG2S = 1998 SS162
                if isdigit(packed[4:6]) and packed[4:6] != "00":
                    desig = (
                        self.__hex2dec(packed[0])
                        + packed[1:3]
                        + " "
                        + packed[3]
                        + packed[-1]
                        + packed[4:6].lstrip("0")
                    )

                if isdigit(packed[4:6]) and packed[4:6] == "00":
                    desig = (
                        self.__hex2dec(packed[0])
                        + packed[1:3]
                        + " "
                        + packed[3]
                        + packed[-1]
                    )

                if not isdigit(packed[4:6]):
                    desig = (
                        self.__hex2dec(packed[0])
                        + packed[1:3]
                        + " "
                        + packed[3]
                        + packed[-1]
                        + self.__hex2dec(packed[4])
                        + packed[5]
                    )

            elif packed[2] == "S":  # ex: T1S3138 = 3138 T-1
                desig = packed[3:] + " " + packed[0] + "-" + packed[1]
        elif packed[0] == "~":
            base62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
            i4 = base62.index(packed[1])
            i3 = base62.index(packed[2])
            i2 = base62.index(packed[3])
            i1 = base62.index(packed[4])
            n = i4 * 62 ^ 3 + i3 * 62 ^ 2 + i2 * 62 ^ 1 + i1 * 62 ^ 0
            desig = f"({n + 620000})"
            # print(f'{packed} {i4} {i3} {i2} {i1} {n} {desig}')
        else:
            print("fail to expand packed designation")

        return desig

    def date_from_designation(self, name: str) -> datetime.datetime:
        isdigit = str.isdigit
        halfmonth_letter = {
            "A": (1, 1),
            "B": (1, 15),
            "C": (2, 1),
            "D": (2, 15),
            "E": (3, 1),
            "F": (3, 15),
            "G": (4, 1),
            "H": (4, 15),
            "J": (5, 1),
            "K": (5, 15),
            "L": (6, 1),
            "M": (6, 15),
            "N": (7, 1),
            "O": (7, 15),
            "P": (8, 1),
            "Q": (8, 15),
            "R": (9, 1),
            "S": (9, 15),
            "T": (10, 1),
            "U": (10, 15),
            "V": (11, 1),
            "W": (11, 15),
            "X": (12, 1),
            "Y": (12, 15),
        }
        if isdigit(name[0:4]) and name[4] == " ":  # start with year
            year = int(name[0:4])
            halfmonth = name[5]
            month, day = halfmonth_letter[halfmonth]
            date = datetime.datetime.strptime(f"{year}-{month}-{day}", "%Y-%m-%d")
        else:
            date = np.nan
        return date

    def date_from_packed_designation(self, packed: str) -> datetime.datetime:
        isdigit = str.isdigit
        if (
            not isdigit(packed[0])
            and isdigit(packed[1:3])
            and (packed[0] in ["I", "J", "K"])
            and len(packed.strip()) == 7
        ):
            try:
                packdt = str(packed).strip()
            except ValueError:
                print("ValueError: Input is not convertable to string.")

            year = self.__hex2dec(packdt[0]) + packdt[1:3]
            halfmonth = float(self.__hex2dec(packdt[3])) - 9
            if halfmonth > 9:
                halfmonth -= 1
            month = str(int(np.ceil(halfmonth / 2)))
            if (halfmonth % 2) == 1:
                day = "01"
            else:
                day = "15"
            result = datetime.datetime.strptime(f"{year}-{month}-{day}", "%Y-%m-%d")
        else:
            result = None
        return result

    def orbit_type(self, a: float, e: float, i: float) -> str:
        """
        Classify asteroid orbit type
        following http://en.wikipedia.org/wiki/Near-Earth_object
        """

        Qt = 1.017
        qt = 0.983
        at = 1
        neo = 1.3
        Q = a * (1 + e)
        q = a * (1 - e)
        t = []

        if q <= neo:
            t.append("NEO")

        if a <= at:
            if Q > qt:
                t.append("Athen")
            else:
                t.append("Atira")
        else:
            if q < Qt:
                t.append("Apollo")
            # Amors (1.0167 < q < 1.3 AU)
            elif Qt < q < neo:
                t.append("Amor")

        # Mars crossers (1.3 < q < 1.6660 AU)
        if neo < q < 1.6660:
            t.append("MarsCrosser")

        # HUNGARIAN Semi-major axis between 1.78 and 2.00 AU. Orbital period of approximately 2.5 years.
        # Low eccentricity of below 0.18. An inclination of 16° to 34°
        if 1.78 <= a <= 2.0 and e <= 0.18 and 16 <= i <= 34:
            t.append("Hungaria")

        # MB:Zona I (2,06-2,5 UA), Zona II (2,5-2,82 UA) y Zona III (2,82-3,28 UA).
        if 2.06 <= a <= 2.5:
            # main belter I
            t.append("MB I")

        if 2.5 <= a <= 2.82:
            # main belter II
            t.append("MB II")

        if 2.82 <= a <= 3.28:
            # main belter II
            t.append("MB III")

        # HILDA: semi-major axis between 3.7 AU and 4.2 AU, an eccentricity less than 0.3, and an inclination less than 20°
        if 3.7 <= a <= 4.2 and e <= 0.3 and i <= 20:
            t.append("Hilda")

        # TNOs 30,103
        if a >= 30.103:
            t.append("TNO")

        # print a,Q,q,t
        t_ = ";".join(t)
        return t_

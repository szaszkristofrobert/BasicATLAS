import re
from pathlib import Path

import numpy as np


class Atlas7V:
    """Python state container and mode-2 reader for ATLAS7V model decks.

    The reader supports the model/control cards written into BasicATLAS
    ``output_synthe.out`` files, including ``READ DECK6`` layer data. It stops
    at ``BEGIN``, matching the mode-2 Fortran reader's return point.
    """

    _ARRAY_SHAPES = {
        "ABROSS": (99,), "TAUROS": (99,),
        "BAL1": (99, 9), "AAL1": (99,), "SAL1": (99,), "XNFPAL": (99, 2), "BAL2": (99, 1),
        "BB1": (99, 7), "XNFPB": (99, 1),
        "BC1": (99, 14), "AC1": (99,), "SC1": (99,), "XNFPC": (99, 2), "BC2": (99, 6),
        "BCA1": (99, 8), "BCA2": (99, 5), "XNFPCA": (99, 2),
        "BFE1": (99, 15), "AFE1": (99,), "SFE1": (99,), "XNFPFE": (99, 1),
        "BHE1": (99, 29), "AHE1": (99,), "SHE1": (99,), "BHE2": (99, 6), "AHE2": (99,),
        "SHE2": (99,), "AHEMIN": (99,), "SIGHE": (99,), "XNFPHE": (99, 3), "XNFHE": (99, 2),
        "BHYD": (99, 8), "AHYD": (99,), "SHYD": (99,), "AH2P": (99,), "BMIN": (99,),
        "AHMIN": (99,), "SHMIN": (99,), "SIGH": (99,), "SIGH2": (99,), "AHLINE": (99,),
        "SHLINE": (99,), "XNFPH": (99, 2), "XNFH": (99,),
        "BK1": (99, 8), "XNFPK": (99, 1),
        "BMG1": (99, 11), "AMG1": (99,), "SMG1": (99,), "XNFPMG": (99, 2), "BMG2": (99, 6),
        "BNA1": (99, 8), "XNFPNA": (99, 1),
        "BO1": (99, 13), "XNFPO": (99, 1), "BO2": (99, 4),
        "BSI1": (99, 11), "ASI1": (99,), "SSI1": (99,), "XNFPSI": (99, 2), "BSI2": (99, 10),
        "DLTDLP": (99,), "HEATCP": (99,), "DLRDLT": (99,), "VELSND": (99,), "GRDADB": (99,),
        "HSCALE": (99,), "FLXCNV": (99,), "VCONV": (99,),
        "ABUND": (99,), "ATMASS": (99,), "XABUND": (99,),
        "FLXERR": (99,), "FLXDRV": (99,), "FLXRAD": (99,),
        "FRESET": (500,), "RCOSET": (500,), "HEIGHT": (99,),
        "IFOP": (20,), "IFPRNT": (15,), "IFPNCH": (15,),
        "ANGLE": (20,), "SURFI": (20,), "KNU": (99,), "PRADK": (99,), "EDENS": (99,),
        "ACCRAD": (99,), "PRAD": (99,), "RHOX": (99,),
        "P": (99,), "XNE": (99,), "XNATOM": (99,), "RHO": (99,), "PTOTAL": (99,),
        "T": (99,), "TKEV": (99,), "TK": (99,), "HKT": (99,), "TLOG": (99,), "HCKT": (99,),
        "VTURB": (99,), "PTURB": (99,),
    }
    _CHAR_ARRAYS = {"ELEM": ((99,), 2), "WORD": ((6,), 1), "FREQID": ((6,), 1), "TITLE": ((74,), 1)}
    _INTEGER_ARRAYS = {"IFOP", "IFPRNT", "IFPNCH"}
    _INTEGER_SCALARS = {
        "NUMCOL", "LETCOL", "LAST", "MORE", "IFFAIL", "MAXPOW", "NULO", "NUHI", "NUMNU",
        "IFWAVE", "IFCORR", "IFPRES", "IFSURF", "IFSCAT", "IFMOL", "NLTEON", "ITER", "NUMITS",
        "NMU", "NRHOX", "ITEMP", "IFCONV", "IFTURB",
    }
    _FLOAT_SCALARS = {
        "MIXLTH", "WTMOLE", "FLUX", "WBEGIN", "DELTAW", "XSCALE", "PZERO", "PCON", "PRADK0",
        "PTURB0", "TEFF", "GRAV", "GLOG", "TRBFDG", "TRBCON", "TRBPOW", "TRBSND",
    }
    _NUMBER = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?")

    def __init__(self):
        for name, shape in self._ARRAY_SHAPES.items():
            dtype = np.int64 if name in self._INTEGER_ARRAYS else np.float64
            setattr(self, name, np.zeros(shape, dtype=dtype))
        for name, (shape, width) in self._CHAR_ARRAYS.items():
            setattr(self, name, np.full(shape, " ", dtype="U{}".format(width)))
        for name in self._INTEGER_SCALARS:
            setattr(self, name, 0)
        for name in self._FLOAT_SCALARS:
            setattr(self, name, 0.0)
        self.WLTE = "LTE "
        self._has_deck6 = False
        self._load_element_data()

    def _load_element_data(self):
        data_file = Path(__file__).resolve().parent / "data" / "solar.csv"
        atomic_numbers = np.loadtxt(data_file, delimiter=",", usecols=(0,))
        masses = np.loadtxt(data_file, delimiter=",", usecols=(1,))
        symbols = np.genfromtxt(data_file, delimiter=",", usecols=(5,), dtype=str)
        for atomic_number, mass, symbol in zip(atomic_numbers, masses, symbols):
            index = int(atomic_number) - 1
            if 0 <= index < 99:
                self.ATMASS[index] = mass
                self.ELEM[index] = symbol

    @classmethod
    def _numbers(cls, text):
        return [float(value.replace("D", "E").replace("d", "e")) for value in cls._NUMBER.findall(text)]

    @staticmethod
    def _set_values(array, values):
        array[:] = 0
        count = min(len(array), len(values))
        array[:count] = values[:count]

    def readin(self, filename, mode=2):
        """Read a supported ATLAS model deck and return this state object.

        ``mode=2`` reads one previously calculated model and returns at its
        ``BEGIN`` record. Mode 1 and streaming mode 20 are not implemented.
        """
        if mode != 2:
            raise ValueError("Only READIN mode=2 is implemented")

        lines = Path(filename).read_text().splitlines()
        line_index = 0
        while line_index < len(lines):
            line = lines[line_index].strip()
            line_index += 1
            if not line:
                continue
            upper = line.upper()
            words = upper.split()
            command = words[0]

            if command == "BEGIN":
                self._finish_mode2()
                return self
            if command == "END":
                raise ValueError("READIN mode=2 does not accept an END card")
            if command == "READ":
                if len(words) > 1 and words[1] in ("DECK", "DECK6"):
                    line_index = self._read_deck6(lines, line_index, line)
                elif len(words) > 1 and words[1] == "MOLECULES":
                    continue
                else:
                    raise ValueError("Unsupported READ card in mode=2: {}".format(line))
            elif command == "TEFF":
                self._read_teff_gravity(line)
            elif command == "TITLE":
                self.TITLE[:] = list(line[6:80].ljust(74)[:74])
            elif command == "OPACITY":
                match = re.search(r"\bIFOP\b(.*)$", line, re.IGNORECASE)
                if match:
                    self._set_values(self.IFOP, [int(value) for value in self._numbers(match.group(1))])
            elif command == "ITERATIONS":
                values = self._numbers(line[len(words[0]):])
                if values:
                    self.NUMITS = int(values[0])
                    self.IFPNCH[:] = 0
                    if self.NUMITS > 0:
                        self.IFPNCH[min(self.NUMITS, len(self.IFPNCH)) - 1] = 1
                self._read_print_punch(line)
            elif command in ("PRINT", "PUNCH"):
                self._read_print_punch(line)
            elif command == "ABUNDANCE":
                self._read_abundances(line)
            elif command == "MOLECULES":
                if len(words) > 1:
                    self.IFMOL = int(words[1] == "ON")
            elif command == "SURFACE":
                self._read_surface(line)
            elif command == "CORRECTION":
                self.IFCORR = int(len(words) > 1 and words[1] == "ON")
                if not self.IFCORR:
                    self.FLXERR[:] = 0
                    self.FLXDRV[:] = 0
            elif command == "PRESSURE":
                self.IFPRES = int(len(words) > 1 and words[1] == "ON")
            elif command == "SCATTERING":
                self.IFSCAT = int(len(words) > 1 and words[1] == "ON")
            elif command == "CONVECTION":
                self._read_convection(line)
            elif command == "TURBULENCE":
                self._read_turbulence(line)
            elif command == "PRADK":
                values = self._numbers(line[len(words[0]):])
                if values:
                    self.PRADK0 = values[0]
            elif command == "LTE":
                self.NLTEON = 0
                self.WLTE = "LTE "
            elif command == "NLTE":
                self.NLTEON = 1
                self.WLTE = "NLTE"
            else:
                raise ValueError("Unsupported mode=2 input card: {}".format(line))

        raise ValueError("Mode=2 input ended before a BEGIN card")

    def _read_teff_gravity(self, line):
        teff = re.search(r"\bTEFF\s+(" + self._NUMBER.pattern + r")", line, re.IGNORECASE)
        gravity = re.search(r"\bGRAVITY\s+(" + self._NUMBER.pattern + r")", line, re.IGNORECASE)
        if teff:
            self.TEFF = float(teff.group(1).replace("D", "E").replace("d", "e"))
            self.FLUX = 5.6697e-5 / 12.5664 * self.TEFF ** 4
        if gravity:
            self.GRAV = float(gravity.group(1).replace("D", "E").replace("d", "e"))
            if self.GRAV < 10.0:
                self.GRAV = 10.0 ** self.GRAV
            self.GLOG = np.log10(self.GRAV)

    def _read_print_punch(self, line):
        for keyword, array in (("PRINT", self.IFPRNT), ("PUNCH", self.IFPNCH)):
            match = re.search(r"\b{}\b(.*?)(?=\b(?:PRINT|PUNCH)\b|$)".format(keyword), line, re.IGNORECASE)
            if match:
                self._set_values(array, [int(value) for value in self._numbers(match.group(1))])

    def _read_abundances(self, line):
        scale = re.search(r"\bSCALE\s+(" + self._NUMBER.pattern + r")", line, re.IGNORECASE)
        if scale:
            self.XSCALE = float(scale.group(1).replace("D", "E").replace("d", "e"))
        change = re.search(r"\bCHANGE\b(.*)$", line, re.IGNORECASE)
        if change:
            values = self._numbers(change.group(1))
            for index in range(0, len(values) - 1, 2):
                atomic_number = int(values[index])
                abundance = values[index + 1]
                if 1 <= atomic_number <= 99:
                    self.ABUND[atomic_number - 1] = abundance

    def _read_surface(self, line):
        words = line.upper().split()
        if len(words) < 2 or words[1] == "OFF":
            self.IFSURF = 0
            self.NMU = 1
        elif words[1] == "FLUX":
            self.IFSURF = 1
            self.NMU = 1
        elif words[1] == "INTENSITY":
            values = self._numbers(line[len("SURFACE INTENSITY"):])
            self.NMU = int(values[0]) if values else 0
            self.ANGLE[:] = 0
            self.ANGLE[:self.NMU] = values[1:self.NMU + 1]
            self.IFSURF = 2
        else:
            raise ValueError("Unsupported SURFACE card: {}".format(line))

    def _read_convection(self, line):
        words = line.upper().split()
        if len(words) < 2 or words[1] == "OFF":
            self.IFCONV = 0
            self.MIXLTH = 1.0
            for name in ("DLTDLP", "HEATCP", "DLRDLT", "VELSND", "GRDADB", "HSCALE", "FLXCNV", "VCONV"):
                getattr(self, name)[:] = 0
        elif words[1] == "ON":
            values = self._numbers(line[len("CONVECTION ON"):])
            self.IFCONV = 1
            self.MIXLTH = values[0] if values else 1.0
        else:
            raise ValueError("Unsupported CONVECTION card: {}".format(line))

    def _read_turbulence(self, line):
        words = line.upper().split()
        if len(words) < 2 or words[1] == "OFF":
            self.IFTURB = 0
            self.TRBFDG = self.TRBPOW = self.TRBSND = self.TRBCON = 0.0
        elif words[1] == "ON":
            values = self._numbers(line[len("TURBULENCE ON"):])
            self.IFTURB = 1
            values += [0.0] * (4 - len(values))
            self.TRBFDG, self.TRBPOW, self.TRBSND, self.TRBCON = values[:4]
        else:
            raise ValueError("Unsupported TURBULENCE card: {}".format(line))

    def _read_deck6(self, lines, line_index, line):
        words = line.split()
        if len(words) < 3:
            raise ValueError("READ DECK6 is missing the layer count: {}".format(line))
        n_layers = int(float(words[2]))
        if not 1 <= n_layers <= 99:
            raise ValueError("READ DECK6 layer count must be between 1 and 99")
        column_text = " ".join(words[3:])
        columns = [name.strip().upper() for name in column_text.split(",") if name.strip()]
        rows = []
        while len(rows) < n_layers and line_index < len(lines):
            row = lines[line_index].strip()
            line_index += 1
            if not row:
                continue
            values = self._numbers(row)
            if len(values) < 7:
                raise ValueError("Invalid READ DECK6 row {}: expected at least 7 values".format(len(rows) + 1))
            rows.append(values)
        if len(rows) != n_layers:
            raise ValueError("READ DECK6 expected {} rows, found {}".format(n_layers, len(rows)))

        self.NRHOX = n_layers
        self._has_deck6 = words[1].upper() == "DECK6"
        if len(columns) < 7:
            columns = ["RHOX", "T", "P", "XNE", "ABROSS", "ACCRAD", "VTURB"] + columns[7:]
        targets = {
            "RHOX": "RHOX", "T": "T", "TEMP": "T", "P": "P", "XNE": "XNE",
            "ABROSS": "ABROSS", "ACCRAD": "PRAD", "PRAD": "PRAD", "VTURB": "VTURB",
            "FLXRAD": "FLXRAD", "FLXCNV": "FLXCNV", "VCONV": "VCONV", "VELSND": "VELSND",
        }
        for row_index, row in enumerate(rows):
            for column_index, value in enumerate(row):
                if column_index < len(columns) and columns[column_index] in targets:
                    getattr(self, targets[columns[column_index]])[row_index] = value
        if self.RHOX[0] < 0.0:
            self.RHOX[:n_layers] = np.power(10.0, self.RHOX[:n_layers])
        return line_index

    @staticmethod
    def _integrate(x, values, start):
        count = len(x)
        result = np.empty(count, dtype=np.float64)
        result[0] = start
        if count == 1:
            return result
        if count < 4:
            for index in range(count - 1):
                result[index + 1] = result[index] + 0.5 * (values[index] + values[index + 1]) * (x[index + 1] - x[index])
            return result

        a = np.zeros(count)
        b = np.zeros(count)
        c = np.zeros(count)
        b[0] = (values[1] - values[0]) / (x[1] - x[0])
        a[0] = values[0] - x[0] * b[0]
        b[-1] = (values[-1] - values[-2]) / (x[-1] - x[-2])
        a[-1] = values[-1] - x[-1] * b[-1]
        for index in range(1, count - 1):
            previous = index - 1
            slope = (values[index] - values[previous]) / (x[index] - x[previous])
            c[index] = (
                values[index + 1] / ((x[index + 1] - x[index]) * (x[index + 1] - x[previous]))
                + (values[previous] / (x[index + 1] - x[previous]) - values[index] / (x[index + 1] - x[index]))
                / (x[index] - x[previous])
            )
            b[index] = slope - (x[index] + x[previous]) * c[index]
            a[index] = values[previous] - x[previous] * slope + x[index] * x[previous] * c[index]
        c[1] = 0.0
        b[1] = (values[2] - values[1]) / (x[2] - x[1])
        a[1] = values[1] - x[1] * b[1]
        c[2] = 0.0
        b[2] = (values[3] - values[2]) / (x[3] - x[2])
        a[2] = values[2] - x[2] * b[2]
        for index in range(1, count - 1):
            if c[index] == 0.0:
                continue
            following = index + 1
            weight = abs(c[following]) / (abs(c[following]) + abs(c[index]))
            a[index] = a[following] + weight * (a[index] - a[following])
            b[index] = b[following] + weight * (b[index] - b[following])
            c[index] = c[following] + weight * (c[index] - c[following])
        a[-2], b[-2], c[-2] = a[-1], b[-1], c[-1]

        for index in range(count - 1):
            left = x[index]
            right = x[index + 1]
            increment = (
                a[index]
                + b[index] / 2.0 * (right + left)
                + c[index] / 3.0 * ((right + left) * right + left * left)
            ) * (right - left)
            result[index + 1] = result[index] + increment
        return result

    def _finish_mode2(self):
        if self.NRHOX == 0:
            raise ValueError("Mode=2 input reached BEGIN before loading a model structure")
        if self._has_deck6:
            count = self.NRHOX
            self.ACCRAD[:count] = self.PRAD[:count]
            self.PRAD[:count] = self._integrate(
                self.RHOX[:count], self.ACCRAD[:count], self.ACCRAD[0] * self.RHOX[0]
            )
            self.PRADK[:count] = self.PRAD[:count] + self.PRADK0

        if self.ABUND[0] < 0.0:
            self.ABUND[0] = 10.0 ** self.ABUND[0]
        if self.ABUND[1] < 0.0:
            self.ABUND[1] = 10.0 ** self.ABUND[1]
        metal_abundances = self.ABUND[2:]
        positive = metal_abundances > 0.0
        metal_abundances[positive] = np.log10(metal_abundances[positive])
        self.XABUND[:2] = self.ABUND[:2]
        self.XABUND[2:] = np.power(10.0, metal_abundances) * self.XSCALE
        self.WTMOLE = float(np.sum(self.XABUND * self.ATMASS))

        count = self.NRHOX
        temperature = self.T[:count]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            self.TK[:count] = 1.38054e-16 * temperature
            self.HKT[:count] = 6.6256e-27 / self.TK[:count]
            self.HCKT[:count] = self.HKT[:count] * 2.99792458e10
            self.TKEV[:count] = 8.6171e-5 * temperature
            self.TLOG[:count] = np.log(temperature)
            self.XNATOM[:count] = self.P[:count] / self.TK[:count] - self.XNE[:count]
            self.RHO[:count] = self.XNATOM[:count] * self.WTMOLE * 1.660e-24
            self.PTURB[:count] = 0.5 * self.RHO[:count] * self.VTURB[:count] ** 2


atlas7v = Atlas7V

__all__ = ["Atlas7V", "atlas7v"]
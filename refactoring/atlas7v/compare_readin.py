"""Compare the Fortran and Python implementations of ATLAS7V READIN.

The Fortran probe is compiled on demand from the checked-in reference source.
It reads the same deck as :class:`Atlas7V`, writes selected COMMON-block
values to a simple text file, and the Python side compares those values.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from atlas7v import Atlas7V


FIELDS = [
    "TEFF", "GRAV", "GLOG", "NUMITS", "NRHOX", "NUMNU",
    "NULO", "NUHI", "IFPRES", "IFSURF", "IFSCAT",
    "IFMOL", "NLTEON", "XSCALE", "PRADK0",
    "RHOX", "T", "P", "XNE", "ABROSS", "PRAD",
    "VTURB", "ABUND", "IFOP", "IFPRNT", "IFPNCH",
]

INTEGER_FIELDS = {
    "NUMITS", "NRHOX", "NUMNU", "NULO", "NUHI", "IFPRES", "IFSURF",
    "IFSCAT", "IFMOL", "NLTEON", "IFOP", "IFPRNT", "IFPNCH",
}
ARRAY_FIELDS = {
    "RHOX", "T", "P", "XNE", "ABROSS", "PRAD", "VTURB", "ABUND",
    "IFOP", "IFPRNT", "IFPNCH",
}
DEPTH_FIELDS = {"RHOX", "T", "P", "XNE", "ABROSS", "PRAD", "VTURB"}


_FORTRAN_PROBE = r"""
      PROGRAM COMPARE_READIN
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KW=99)
      CHARACTER*256 DECK,MOLECULES,DUMP
      INTEGER I,MODE
      INTEGER NUMITS,NRHOX,NUMNU,NULO,NUHI,IFPRES,IFSURF,IFSCAT
      INTEGER IFMOL,NLTEON,IFOP(20),IFPRNT(15),IFPNCH(15)
      REAL*8 TEFF,GRAV,GLOG,XSCALE,PRADK0
      REAL*8 RHOX(KW),T(KW),P(KW),XNE(KW),ABROSS(KW),PRAD(KW)
      REAL*8 VTURB(KW),ABUND(99)
      COMMON /ABROSS/ABROSS,TAUROS(99)
      COMMON /ELEM/ABUND,ATMASS(99),ELEM(99)
      COMMON /FRESET/FRESET(500),RCOSET(500),NULO,NUHI,NUMNU,IFWAVE,
     1 WBEGIN,DELTAW
      COMMON /IF/IFCORR,IFPRES,IFSURF,IFSCAT,IFMOL,NLTEON,IFOP
      COMMON /ITER/ITER,IFPRNT,IFPNCH,NUMITS
      COMMON /JUNK/XSCALE,TITLE(74),FREQID(6),WLTE
      COMMON /PZERO/PZERO,PCON,PRADK0,PTURB0,KNU(99),PRADK(99),EDENS(99)
      COMMON /RAD/ACCRAD(99),PRAD
      COMMON /RHOX/RHOX,NRHOX
      COMMON /STATE/P,XNE,XNATOM(99),RHO(99),PTOTAL(99)
      COMMON /TEFF/TEFF,GRAV,GLOG
      COMMON /TEMP/T,TKEV(99),TK(99),HKT(99),TLOG(99),HCKT(99),ITEMP
      COMMON /TURBPR/VTURB,PTURB(99),TRBFDG,TRBCON,TRBPOW,TRBSND,IFTURB

      CALL GET_COMMAND_ARGUMENT(1,DECK)
      CALL GET_COMMAND_ARGUMENT(2,MOLECULES)
      CALL GET_COMMAND_ARGUMENT(3,DUMP)
      OPEN(5,FILE=TRIM(DECK),STATUS='OLD',ACTION='READ')
      OPEN(2,FILE=TRIM(MOLECULES),STATUS='OLD',ACTION='READ')
      OPEN(9,FILE=TRIM(DUMP),STATUS='REPLACE',ACTION='WRITE')
      MODE=2
      CALL BLOCKR
      CALL READIN(MODE)
      WRITE(9,'(A,1X,ES24.16)') 'S TEFF',TEFF
      WRITE(9,'(A,1X,ES24.16)') 'S GRAV',GRAV
      WRITE(9,'(A,1X,ES24.16)') 'S GLOG',GLOG
      WRITE(9,'(A,1X,I12)') 'I NUMITS',NUMITS
      WRITE(9,'(A,1X,I12)') 'I NRHOX',NRHOX
      WRITE(9,'(A,1X,I12)') 'I NUMNU',NUMNU
      WRITE(9,'(A,1X,I12)') 'I NULO',NULO
      WRITE(9,'(A,1X,I12)') 'I NUHI',NUHI
      WRITE(9,'(A,1X,I12)') 'I IFPRES',IFPRES
      WRITE(9,'(A,1X,I12)') 'I IFSURF',IFSURF
      WRITE(9,'(A,1X,I12)') 'I IFSCAT',IFSCAT
      WRITE(9,'(A,1X,I12)') 'I IFMOL',IFMOL
      WRITE(9,'(A,1X,I12)') 'I NLTEON',NLTEON
      WRITE(9,'(A,1X,ES24.16)') 'S XSCALE',XSCALE
      WRITE(9,'(A,1X,ES24.16)') 'S PRADK0',PRADK0
      DO I=1,99
         WRITE(9,'(A,1X,I4,1X,ES24.16)') 'A RHOX',I,RHOX(I)
         WRITE(9,'(A,1X,I4,1X,ES24.16)') 'A T',I,T(I)
         WRITE(9,'(A,1X,I4,1X,ES24.16)') 'A P',I,P(I)
         WRITE(9,'(A,1X,I4,1X,ES24.16)') 'A XNE',I,XNE(I)
         WRITE(9,'(A,1X,I4,1X,ES24.16)') 'A ABROSS',I,ABROSS(I)
         WRITE(9,'(A,1X,I4,1X,ES24.16)') 'A PRAD',I,PRAD(I)
         WRITE(9,'(A,1X,I4,1X,ES24.16)') 'A VTURB',I,VTURB(I)
         WRITE(9,'(A,1X,I4,1X,ES24.16)') 'A ABUND',I,ABUND(I)
      END DO
      DO I=1,20
         WRITE(9,'(A,1X,I4,1X,I12)') 'AI IFOP',I,IFOP(I)
      END DO
      DO I=1,15
         WRITE(9,'(A,1X,I4,1X,I12)') 'AI IFPRNT',I,IFPRNT(I)
         WRITE(9,'(A,1X,I4,1X,I12)') 'AI IFPNCH',I,IFPNCH(I)
      END DO
      CLOSE(9)
      END
"""


def _compile_probe(source: Path, executable: Path, compiler: str | None = None) -> None:
    compiler = compiler or os.environ.get("FC")
    if compiler is None:
        compiler = next(
            (candidate for candidate in ("gfortran", "ifx", "ifort") if shutil.which(candidate)),
            None,
        )
    if compiler is None:
        raise RuntimeError(
            "A Fortran compiler is required (set FC or install gfortran/ifx/ifort)"
        )
    compiler_name = Path(compiler).name.lower()
    if compiler_name == "gfortran":
        flags = [
            "-O0", "-std=legacy", "-ffixed-line-length-none",
            "-ffunction-sections", "-fdata-sections", "-Wl,--gc-sections",
        ]
    else:
        flags = ["-O0", "-fixed"]
    command = [
        compiler, *flags, "-o", str(executable), str(source),
        str(Path(__file__).with_name("atlas7v.for")),
    ]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError as error:
        raise RuntimeError(f"Fortran compiler not found: {compiler}") from error
    except subprocess.CalledProcessError as error:
        diagnostic = (error.stdout + "\n" + error.stderr).strip()
        raise RuntimeError("Could not compile the Fortran READIN probe:\n" + diagnostic) from error


def _fortran_state(
    deck: Path, molecules: Path, compiler: str | None = None
) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="atlas7v-readin-") as directory:
        directory_path = Path(directory)
        source = directory_path / "probe.f"
        executable = directory_path / "probe"
        dump = directory_path / "state.txt"
        source.write_text(_FORTRAN_PROBE)
        _compile_probe(source, executable, compiler)
        result = subprocess.run(
            [str(executable), str(deck), str(molecules), str(dump)],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode:
            raise RuntimeError(
                "Fortran READIN failed with exit code {}:\n{}{}".format(
                    result.returncode, result.stdout, result.stderr
                )
            )
        return _parse_dump(dump)


def _parse_dump(filename: Path) -> dict[str, object]:
    state: dict[str, object] = {}
    for line in filename.read_text().splitlines():
        kind, name, index, value = (line.split() + [None, None])[:4]
        if kind == "S":
            state[name] = float(index)
        elif kind == "I":
            state[name] = int(index)
        elif kind == "A":
            state.setdefault(name, np.zeros(99, dtype=np.float64))[int(index) - 1] = float(value)
        elif kind == "AI":
            state.setdefault(name, np.zeros(20 if name == "IFOP" else 15, dtype=np.int64))[int(index) - 1] = int(value)
    return state


def _python_state(deck: Path) -> dict[str, object]:
    atlas = Atlas7V().readin(deck, mode=2)
    return {name: getattr(atlas, name).copy() if name in ARRAY_FIELDS else getattr(atlas, name)
            for name in FIELDS}


def compare(
    deck: Path, molecules: Path, rtol: float, atol: float, compiler: str | None = None
) -> int:
    reference = _fortran_state(deck, molecules, compiler)
    candidate = _python_state(deck)
    mismatches: list[str] = []
    for name in FIELDS:
        expected = reference[name]
        actual = candidate[name]
        if name in INTEGER_FIELDS:
            if not np.array_equal(actual, expected):
                mismatches.append(f"{name}: expected {expected!r}, got {actual!r}")
            continue
        if name in DEPTH_FIELDS:
            count = int(reference["NRHOX"])
            expected = np.asarray(expected)[:count]
            actual = np.asarray(actual)[:count]
        if not np.allclose(actual, expected, rtol=rtol, atol=atol, equal_nan=True):
            difference = np.max(np.abs(np.asarray(actual) - np.asarray(expected)))
            mismatches.append(f"{name}: max absolute difference {difference:.6e}")
    if mismatches:
        print("READIN comparison FAILED")
        print("\n".join("  " + mismatch for mismatch in mismatches))
        return 1
    print(f"READIN comparison PASSED: {deck}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("deck", type=Path, help="ATLAS input deck to parse")
    parser.add_argument(
        "--molecules",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "data" / "atlas_files" / "molecules.dat",
        help="Fortran unit-2 molecule table (default: BasicATLAS data table)",
    )
    parser.add_argument("--rtol", type=float, default=1e-12)
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument(
        "--compiler",
        help="Fortran compiler executable (default: FC, then gfortran, ifx, or ifort)",
    )
    arguments = parser.parse_args()
    return compare(
        arguments.deck, arguments.molecules, arguments.rtol, arguments.atol, arguments.compiler
    )


if __name__ == "__main__":
    raise SystemExit(main())

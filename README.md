# BasicATLAS

This project is a simple Python wrapper for the **ATLAS** LTE stellar modelling code by [Robert L Kurucz](http://kurucz.harvard.edu/) as well as satellite utilities including the spectral synthesis library **SYNTHE** and opacity distribution (ODF) calculator **DFSYNTHE**. Refer to the literature references below:

* [Kurucz (1970)](https://ui.adsabs.harvard.edu/abs/1970SAOSR.309.....K/abstract): detailed description of an older version of the code (version 5)
* [Kurucz et al. (1974)](https://ui.adsabs.harvard.edu/abs/1974bmae.book.....K/abstract): more on opacity distribution functions
* [Kurucz & Avrett (1981)](https://ui.adsabs.harvard.edu/abs/1981SAOSR.391.....K/abstract): standard reference for the **SYNTHE** code
* [Castelli & Kurucz (2003)](https://ui.adsabs.harvard.edu/abs/2003IAUS..210P.A20C/abstract): description of the "new" ODF format, currently used by **ATLAS**
* [Castelli (2005)](https://ui.adsabs.harvard.edu/abs/2005MSAIS...8...34C/abstract): user guide to **DFSYNTHE**
* [Kurucz (2005a)](https://ui.adsabs.harvard.edu/abs/2005MSAIS...8...14K/abstract): more information on modern versions of **ATLAS** as well as **SYNTHE**
* [Kurucz 2005b](https://ui.adsabs.harvard.edu/abs/2005MSAIS...8...86K/abstract): bound-bound (line) opacity treatment in the code and line lists
* [Kurucz (2014)](https://ui.adsabs.harvard.edu/abs/2014dapb.book...39K/abstract): brief overview of the modern versions of the code (9 and 12)

Also see the websites of [Robert L Kurucz](http://kurucz.harvard.edu/programs.html) and [Fiorella Castelli](https://wwwuser.oats.inaf.it/castelli/).

**ATLAS 9** and **ATLAS 12** differ in the opacity sampling method (opacity distribution functions vs direct sampling) which results in the former being considerably faster and, in principle, somewhat less accurate although (to my understanding) the significance of the uncertainty introduced by the chosen opacity treatment is not well established. **BasicATLAS** works with **ATLAS 9** (support for **ATLAS 12** may be added in the future).

If you are using **BasicATLAS** in your research, please cite [Larkin et al. (2022)](https://arxiv.org/abs/2210.09185)

## Contributors

* Roman Gerasimov (University of Notre Dame, University of California San Diego)
* Mikaela Larkin (University of California San Diego)
* Tianxing Zhou (University of California San Diego)
* Philipp Edelmann (Los Alamos National Laboratory) @ [AAS 242](https://aas.org/meetings/aas242) Hack Together
* Paul Barrett (George Washington University) @ [AAS 242](https://aas.org/meetings/aas242) Hack Together
* Efrain Alvarado (University of California San Diego)
* Evan Kirby (University of Notre Dame)

## Installation

### Grid of restart files

**BasicATLAS** uses a grid of pre-computed model atmospheres as initial guesses for temperature-pressure profiles in structure calculations. This repository contains a reduced grid of models; however, we recommend downloading the full grid from [Google Drive](https://drive.google.com/file/d/1xBhLEdUBZTjtHHg110FVH6G-rYFaPHTk/view?usp=sharing) and placing it in `BasicATLAS/restarts/light.h5`.

### Download legacy code and data files

Run the wget-based bash script `download.com` to retrieve the original **ATLAS**/**SYNTHE**/**DFSYNTHE** Fortran codes and required data files (line lists, chemical constants etc):

```bash
source download.com
```

### Compile the code

To compile the legacy code, both Intel and GNU Fortran compilers are required. Intel compilers are typically not available by default, but can be installed as part of the [oneAPI Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html). Follow the instructions on the oneAPI website to install the Toolkit for your operating system (both GUI and command line installers are available). Once installed, the Toolkit must be initialized before every use by running the `setvars.sh` script:

```bash
sh ~/intel/oneapi/setvars.sh
```

Once installed and initialized, you must be able to use the `ifx` command in your terminal. Then run the compile script in **BasicATLAS**:

```bash
source compile.com
```

Finally, verify that the installation completed successfully by running

```bash
python test.py
```

The test is clean if no output is produced.

### PyTLAS sub-module

**BasicATLAS** is shipped with the **PyTLAS** sub-module, which implements a modified version of the **SYNTHE** suite that runs entirely in memory. **PyTLAS** is necessary e.g. for abundance fitting with [chemfit](https://github.com/Roman-UCSD/chemfit). To enable **PyTLAS**:

```
cd PyTLAS
source compile.com
```

This should generate `spectrv.so`, `synthe.so` and `xnfpelsyn.so` in `BasicATLAS/PyTLAS/bin/`.

### Compare the Python and Fortran READIN implementations

The refactoring parser can be checked against the original Fortran `READIN`
routine with the comparison harness. It compiles a small probe from
`refactoring/atlas7v/atlas7v.for`, parses the same deck with
`refactoring/atlas7v/atlas7v.py`, and compares the requested scalar, control,
and atmospheric fields. Integer/control fields use exact equality; floating
point fields use `rtol=1e-12` and `atol=1e-12` by default.

From the repository root, run:

```bash
python refactoring/atlas7v/compare_readin.py \
  /home/bbzze2/ATLAS_sun_model/output_synthe.out
```

The Fortran probe opens `data/atlas_files/molecules.dat` as unit 2 by default.
Use `--molecules` for another table, `--compiler` to select a compiler, or
`--rtol`/`--atol` to change numerical tolerances. The harness requires a
Fortran compiler (`gfortran`, `ifx`, or `ifort`) in `PATH`, or via `FC`.

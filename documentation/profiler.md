# Profiling the SPECTRV executable

## Starting point

This guide starts from the completed ATLAS run in:

```bash
cd /home/bbzze2/ATLAS_sun_model
```

The folder contains the expected ATLAS outputs:

```text
atlas_control.com
atlas_control_end.com
atlas_control_start.com
output_last_iteration.out
output_main.out
output_summary.out
restart.dat
```

The `ATLAS` stage has finished, but the SYNTHE stage has not. The next step is to generate `output_synthe.out` and its launch script so that `SPECTRV` can read the model.

See `BasicATLAS/atlas.py` for the Python implementation that creates these SYNTHE inputs from the completed ATLAS outputs.

## Purpose

This guide explains how to build and run a separate, instrumented `SPECTRV` executable for function-level profiling with `gprof`.

The profiling build is kept separate from the normal production executable so that:

- the production `bin/spectrv.exe` remains unchanged;
- the instrumented build can use profiling flags without affecting normal runs;
- profiling artifacts remain under `BasicATLAS/bin/profile/`; and
- each synthesis batch can retain its own `gmon.out` file.

The profiler targets `SPECTRV` because it is the final computational stage of the SYNTHE pipeline, reads the completed atmospheric model, and performs the expensive line-by-line spectrum calculation.

## Requirements

Install or make available:

- GNU Fortran (`gfortran`)
- GNU binutils, including `gprof`
- Python 3, for the BasicATLAS Python wrapper
- The existing BasicATLAS source, data files, and model inputs

On Linux, verify the commands before building:

```bash
command -v gfortran
command -v gprof
command -v python3
```

## 1. Build the instrumented executable

From the `BasicATLAS` directory, build the profiled executable without modifying the normal executable:

```bash
mkdir -p bin/profile

gfortran -fno-automatic -w -O3 -pg -c src/spectrv.patched.for -std=legacy

gfortran -pg -o bin/profile/spectrv.profile.exe spectrv.patched.o atlas7v.o -std=legacy

rm spectrv.patched.o atlas7v.o
```

The `-pg` option causes the compiler to insert profiling calls and makes the generated program produce `gmon.out` when it exits normally or via `EXIT`.

`-Winvalid-pch`, `-fstack-protector`, or other additional flags are not required for this workflow. Keep the same compile settings as the production build as much as possible so that the profile reflects the normal optimized code.

If `atlas7v.o` is not already available, build it first:

```bash
mkdir -p bin/profile

gfortran -fno-automatic -w -O3 -c src/atlas7v.for -std=legacy
```

The current repository build also produces the unprofiled executable with:

```bash
gfortran -fno-automatic -w -O3 -c src/spectrv.patched.for -std=legacy
gfortran spectrv.patched.o atlas7v.o -o bin/spectrv.exe -std=legacy
rm atlas7v.o spectrv.patched.o
```

## 2. Generate the SYNTHE inputs

The ATLAS outputs are already available in the current directory. Generate the `SPECTRV` input and launch script from them with the Python wrapper:

```bash
cd /home/bbzze2/ATLAS_sun_model
python3 /home/bbzze2/BasicATLAS/atlas.py synthe --output-dir /home/bbzze2/ATLAS_sun_model --min-wl 2850 --max-wl 3000 --res 600000 --vturb 1.5
```

The command creates:

```text
output_synthe.out
synthe_launch.com
synthe_1/
synthe_2/
```

The command must be run before trying to profile `SPECTRV`. If the exact Python CLI is unavailable, use the project’s existing profile notebook or equivalent Python call that invokes `atlas.synthe(...)` with the same run directory and wavelength range.

The generated `output_synthe.out` is the input that `SPECTRV` reads through standard input. `synthe_launch.com` is the full SYNTHE launch script and currently contains the command used for each executable.

## 3. Replace the SPECTRV command in the launcher

The profiling launch must use the instrumented binary:

```bash
/home/bbzze2/BasicATLAS/bin/profile/spectrv.profile.exe
```

Change the `SPECTRV` command inside `synthe_launch.com` from the normal executable to the profiled executable:

```bash
sed -i 's|/home/bbzze2/BasicATLAS/bin/spectrv.exe|/home/bbzze2/BasicATLAS/bin/profile/spectrv.profile.exe|' \
  /home/bbzze2/ATLAS_sun_model/synthe_launch.com
```

Verify the edited command:

```bash
grep -n "spectrv" /home/bbzze2/ATLAS_sun_model/synthe_launch.com
```

Do not use the command that redirects `output_synthe.out` directly unless you also provide the required `fort.25` parameters. The full launched SYNTHE script is safer for a normal profile run because it creates the required `fort.*` files and launches `SPECTRV` with its expected input.

## 4. Run the profiled synthesis

Run one batch from its own working directory so that its profiling data is not overwritten by another batch:

```bash
cd /home/bbzze2/ATLAS_sun_model/synthe_1
```

Then run the generated launcher:

```bash
bash ../synthe_launch.com
```

The launcher creates the required `fort.*` files and runs the entire SYNTHE/SPECTRV batch. It writes `gmon.out` to the batch directory when the profiled `SPECTRV` command finishes.

The direct command is only appropriate when the single `SPECTRV` stage is being tested. It still requires `fort.25` and the appropriate data files:

```bash
cd /home/bbzze2/ATLAS_sun_model/synthe_1
ln -sf /home/bbzze2/BasicATLAS/data/synthe_files/molecules.dat fort.2
cat <<'EOF' > fort.25
0.0       0.        1.        0.        0.        0.        0.        0.
0.
RHOXJ     R1        R101      PH1       PC1       PSI1      PRDDOP    PRDPOW
EOF

/home/bbzze2/BasicATLAS/bin/profile/spectrv.profile.exe \
  < /home/bbzze2/ATLAS_sun_model/output_synthe.out \
  > spectrv.out
```

The `fort.25` file is not part of the completed ATLAS model; it supplies the eight `SPECTRV` control parameters at the start of the program.

## 5. Verify that profiling ran

After the executable exits, confirm that `gmon.out` exists:

```bash
ls -lh gmon.out
file gmon.out
```

The file must be non-empty. If `gmon.out` is absent, check that:

- the executable was built with `-pg`;
- the program terminated normally;
- the executable is running from the expected directory;
- the profile target is the profiled executable; and
- no earlier process or shell command deleted the output file.

## 6. Generate the `gprof` report

Run `gprof` from the directory containing the executable and `gmon.out`:

```bash
cd /home/bbzze2/ATLAS_sun_model/synthe_1

gprof /home/bbzze2/BasicATLAS/bin/profile/spectrv.profile.exe \
  gmon.out > gprof_spectrv.txt
```

The report is written to `gprof_spectrv.txt`. Use the same directory for the executable and the data file only if the executable is installed there; otherwise, provide the full path to the executable.

To produce a separate report for each synthesis batch:

```bash
for batch in synthe_1 synthe_2; do
  cd "/home/bbzze2/ATLAS_sun_model/$batch"
  gprof /home/bbzze2/BasicATLAS/bin/profile/spectrv.profile.exe \
    gmon.out > gprof_spectrv.txt
  cd ..
done
```

## 7. Interpret the report

Open the generated report:

```bash
less gprof_spectrv.txt
```

The most useful sections are:

- **Flat profile**: cumulative time, self time, and percentage for each routine;
- **Call graph**: callers, callees, and call counts;
- **Index**: the routine names used by the call graph and flat profile.

Focus first on routines with:

1. a high percentage of total self time;
2. many calls to expensive child routines; and
3. poor locality or repeated work in a loop.

For example, look for entries such as `JOSH`, `READIN`, `READMOL`, line-opacity routines, or numerical interpolation routines. The precise routines depend on the selected `SPECTRV` path and the model configuration.

A routine with high total time but low self time may be mainly a caller of expensive work. A routine with high self time is the more likely target for direct optimization.

Use call counts and function names together. A routine can appear expensive because it is called many times, even when each call is inexpensive.

## 8. Keep profiling artifacts isolated

Store all profile outputs in the project's profile directory:

```text
BasicATLAS/
├── bin/
│   └── profile/
│       ├── spectrv.profile.exe
│       ├── spectrv.profile.o
│       └── gmon.out
└── documentation/
    └── profiler.md
```

When running multiple model batches, save each batch's `gmon.out` or `gprof` report separately, for example:

```text
ATLAS_sun_model/
├── synthe_1/
│   ├── gmon.out
│   ├── gprof_spectrv.txt
│   └── spectrv.out
└── synthe_2/
    ├── gmon.out
    ├── gprof_spectrv.txt
    └── spectrv.out
```

Remove old profiling outputs before a new measurement if you need a clean run:

```bash
rm -f gmon.out gprof_spectrv.txt spectrv.out
```

## 9. Troubleshooting

### `gmon.out` is missing

Rebuild with `-pg` and verify that the executable is the profiler-targeted binary:

```bash
file bin/profile/spectrv.profile.exe
strings bin/profile/spectrv.profile.exe | grep gmon
```

The executable should contain profiling instrumentation, and it should file-check as an executable.

### `gprof` reports `No profile data`

Check the current working directory and verify that `gmon.out` was generated there:

```bash
pwd
ls -lh gmon.out
```

`gmon.out` must be generated in the same directory from which `gprof` is run, or its absolute path must be provided to `gprof`.

### `gprof` cannot find the executable

Use the executable's full path:

```bash
gprof /home/bbzze2/BasicATLAS/bin/profile/spectrv.profile.exe gmon.out
```

### The executable exits immediately

Check `spectrv.out` and the standard error stream for the underlying Fortran error. The profiled executable must still receive the complete `SPECTRV` input and all required `fort.*` files.

### The normal executable is used

Verify the command in the launcher:

```bash
grep -n "spectrv" /home/bbzze2/ATLAS_sun_model/synthe_launch.com
```

The launcher should reference `BasicATLAS/bin/profile/spectrv.profile.exe` for profiling runs.

## Recommended workflow

1. Start from the completed ATLAS run in `ATLAS_sun_model`.
2. Generate the missing SYNTHE inputs with `atlas.synthe(...)`.
3. Replace the `SPECTRV` command in `synthe_launch.com` with the profiled executable.
4. Run the launcher from `synthe_1` or another dedicated batch directory.
5. Confirm that `gmon.out` exists.
6. Generate a `gprof` report.
7. Identify high-self-time routines first.
8. Apply one focused optimization or refactoring.
9. Rebuild and compare the new profile with the baseline.

This sequence starts from the current model state and keeps the ATLAS/SYNTHE stages separated.

# Release notes

## 2.2.0

### Added

- Windows wheels, for CPython 3.10 through 3.14 on AMD64. These are the first
  Windows wheels in several years.
- aarch64 Linux wheels, both manylinux and musllinux.
- API documentation published at <https://melizalab.github.io/libtfr/>, with the
  python API rendered by pdoc and the C API by doxygen (#11).
- `SPDX-License-Identifier` headers throughout the source.

### Changed

- **C API.** `cmplx_t` is now `double[2]` rather than `double _Complex`, and
  `mtcomplex()` and `mtm_zspec()` take `cmplx_t *`. The binary layout is
  identical, so the ABI is unchanged and compiled code keeps working, but C
  sources that perform native complex arithmetic on these values will no longer
  compile and need to index `[0]` and `[1]` instead. This is what makes MSVC,
  which does not implement C99 `_Complex`, able to build the library.
- **Declared license corrected to `GPL-2.0-or-later`.** The package metadata had
  said `GPL-3.0-or-later` while COPYING and every source header said version 2.
  The terms of the license have not changed; only the declaration was wrong.
  `pyproject.toml` now uses the PEP 639 `license` and `license-files` fields.
- Linux wheels build against `manylinux_2_28`.
- `FFTW_NO_Complex` is defined on every platform rather than only on Windows, so
  there is a single code path for the complex arithmetic instead of a
  Windows-only one that CI exercised less often.
- Python docstrings converted to Google style, and several corrected: `dpss`
  documented parameters it does not have, `log_fgrid` documented a `base`
  argument that does not exist, and `mtstft_pt` documented one return value
  where it returns two.

### Removed

- The MATLAB/mex interface (`matlab/`), which was untested and unreferenced.
- The vendored scons tools in `site_scons/`. Doxygen is now invoked directly
  from `Sconstruct`, and the python extension has been built by setuptools
  rather than scons for a long time.
- `test/profile.py`, whose profiling scaffolding was commented out and could not
  have worked against the current API.

### Fixed

- A numpy 2.5 deprecation in `mfft_precalc`, which assigned to `.shape` in place
  and would have become an error in a later release.
- `mtcomplex()` and `mtm_zspec()` were defined with signatures that did not
  match their declarations in `tfr.h`.
- `scons docs` produced nothing: the doxygen configuration pointed at `tfr.h`
  rather than `include/tfr.h`, so it had no input.
- `LIBTFR_VERSION` in `tfr.h` had been stuck at `2.0.0`.
- `examples/tfr_tm.py` had not run since the module-level transform functions
  became methods on `mfft`. It now uses the current API and checks the
  reassigned spectrogram against the known frequency law of its test signal.
- The C test program ran its transforms over uninitialized memory, and exited 0
  unconditionally, so the C library CI job could not fail. It now asserts and
  covers `mtcomplex()` and `mtm_zspec()`, which nothing else tests.
- The published documentation site is no longer advertised as unsupported on
  Windows, and the sdist ships the files `Sconstruct` needs.

---

Releases before 2.2.0 predate this file; see the git history for those.

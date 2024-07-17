# CelluloseBuilder

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cgtbatista.github.io/CelluloseBuilder.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cgtbatista.github.io/CelluloseBuilder.jl/dev/)
[![Build Status](https://github.com/cgtbatista/CelluloseBuilder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cgtbatista/CelluloseBuilder.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cgtbatista/CelluloseBuilder.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cgtbatista/CelluloseBuilder.jl)
[![Coverage](https://coveralls.io/repos/github/cgtbatista/CelluloseBuilder.jl/badge.svg?branch=main)](https://coveralls.io/github/cgtbatista/CelluloseBuilder.jl?branch=main)

A Julia package that reimplements the Cellulose Builder software, aiming to simplify the code and provide more flexibility for dealing with different crystallographic parameters or even newest microfibril models. The package generates cartesian coordinates for a given crystalographic data of the specified structure in XYZ and PDB formats, suitable for using as starting configurations in molecular dynamics simulations and other calculations. The crystaline polymorphs Iα, Iβ, II, and III(I) of practically any size can be constructed, including parallelepipeds, plant cell wall cellulose elementary fibrils of any length, and monolayers. Periodic boundary conditions (PBC) along the crystallographic directions are easily imposed and it can consider covalent bounding through the PBC. The package also generates atom connectivity files in PSF format, required by well-known simulation packages such as NAMD, CHARMM, and others. CelluloseBuilder.jl is written in Julia and should run on any platform that supports Julia and both are freely available.


## Original Script Update

The original [Cellulose Builder script](https://code.google.com/archive/p/cellulose-builder/) is written in bash, Perl, and Tcl, requiring additional software like Octave to run on a personal computer (Gomes and Skaf, 2012). There is an [online version](http://cces-sw.iqm.unicamp.br/cces/admin/cellulose/view;jsessionid=) that is easy to use, however it does not allow adaptations on the source script. Indeed, this julia package still requires the VMD software (psfgen) to work (Humphrey, 1996). Beyond the PSF, other topology formats will be implemented later using TopoTools.

This package, CelluloseBuilder.jl, aims to:

* Simplify the code
* Provide more flexibility for dealing with different parameters
* Implement new features not available in the original version

## Documentation

For more information, please check out our [stable documentation](https://cgtbatista.github.io/CelluloseBuilder.jl/stable/) or [dev documentation](https://cgtbatista.github.io/CelluloseBuilder.jl/dev/).

## Build Status and Coverage

Check out our [build status](https://github.com/cgtbatista/CelluloseBuilder.jl/actions/workflows/CI.yml?query=branch%3Amain) and [code coverage](https://codecov.io/gh/cgtbatista/CelluloseBuilder.jl) reports.

## References

T. C. F. Gomes, M. S. Skaf, J. Comput. Chem. 2012, 33, 1338-1346. [10.1002/jcc.22959](https://onlinelibrary.wiley.com/doi/10.1002/jcc.22959).

W. Humphrey, A. Dalke, K. Schulten, J. Mol. Graph. 1996, 14, 33-38. [10.1016/0263-7855(96)00018-5](https://linkinghub.elsevier.com/retrieve/pii/0263785596000185).

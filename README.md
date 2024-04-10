[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3516436.svg)](https://doi.org/10.5281/zenodo.3516436)

![missing TInF logo][logo]

[logo]: https://nheri-simcenter.github.io/TurbulenceInflowTool/docs/NHERI-TInF-icon.png "Turbulence Inflow Tool Logo"

# TurbulenceInflowTool

The Turbulence Inflow Tool (TInF) is designed to collect all required properties and parame- ters needed for various turbulence inflow models in OpenFOAM, and to augment an existing wind-around-a-building model by adding the necessary sections to respective parameter definition files.

The generic workflow involved is as follows.
1. Build your OpenFOAM model as you would without using a turbulence inflow model. Use a generic patch with a suitable name for you will need to identify that patch by its name inside TInF.
2. Run TInF following, identify your model folder, adjust the parameters as desired, and export to your model definition. Consult Chapter 4 for details on those steps.
3. Run OpenFOAM using the updated model definition.
The tool also provides a Save to file and Open from file functionality that will allow you to define and share complex sets of settings and parameters for the supported turbulence inflow models and, such, efficiently and reliably apply the same parameters to several different models.

Technical detail is available from the TInF manual available through: [https://nheri-simcenter.github.io/TinF-Documentation/](https://nheri-simcenter.github.io/TinF-Documentation/)


# references

1. [https://nheri-simcenter.github.io/TurbulenceInflowTool/](https://nheri-simcenter.github.io/TurbulenceInflowTool/)
2. [https://simcenter.designsafe-ci.org/research-tools/tinf/](https://simcenter.designsafe-ci.org/research-tools/tinf/)

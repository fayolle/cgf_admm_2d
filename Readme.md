# Iterative distance approximation 
Implementation in Matlab (and in 2D) of the method described in Section 5 of this [paper](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.12611). Two variants are provided: One implementation is using finite differences, and is suitable for data defined on a regular grid, such as images (Euclidean Distance Transform, skeleton computation). The second implementation is using FEM and is suitable for computing the distance to a boundary (curve) on a triangulated data. 


## Examples 
- 'demo' in the directory 'cgf_admm_grid' for a finite difference based approach. It uses Matlab's builtins functions for visualizing the results. 
- 'demo' in the directory 'cgf_admm_FEM' for a FEM based approach. It uses Matlab's builtins functions for visualizing the results and generates VTK files as well (legacy, text based format) that can be visualized with [paraview](https://www.paraview.org/). 


## Input data
- 'demo' in the directory 'cgf_admm_grid' expects black and white images. Example images are provided in the directory 'data'. 
- 'demo' in the directory 'cgf_admm_FEM' expects as input a triangulated 2D domain. The file format exported by the program [triangle](https://www.cs.cmu.edu/~quake/triangle.html) is used. Examples of data for two simple domains (a square and a disk), and a more complicated one (a rider) are provided in the directory 'data'. 


## Reference 
Link to the [paper](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.12611) where the method was introduced. The corresponding bibtex entry is 
```
@article{Belyaev2015,
author = {Belyaev, Alexander G. and Fayolle, Pierre-Alain},
title = {On Variational and PDE-Based Distance Function Approximations},
journal = {Computer Graphics Forum},
volume = {34},
number = {8},
pages = {104-118},
year = {2015}
}

```

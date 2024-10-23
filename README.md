# Leaky Hockey Stick Loss: The First Negatively Divergent Margin-based Loss Function for Classification

This repository provides an R package, `lhsc`, for fitting leaky hockey stick classifier, from the paper: 

Oh-Ran Kwon and Hui Zou (2023) **Leaky Hockey Stick Loss: The First Negatively Divergent Margin-based Loss Function for Classification**, *Journal of Machine Learning Research*, 24(239), 1-40. URL: [http://jmlr.org/papers/v24/22-1104.html](http://jmlr.org/papers/v24/22-1104.html)

## Installation of the R package `lhsc`
You may need the Fortran compiler (gfortran) on your system, as this R package contains Fortran code. To proceed, open an `R` prompt and type the following commands:
```
library(devtools)
install_github("ohrankwon/lhsc")
library(lhsc)
```

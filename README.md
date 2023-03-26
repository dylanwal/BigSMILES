# BigSMILES Parser

---
---

![PyPI](https://img.shields.io/pypi/v/bigsmiles)
![downloads](https://static.pepy.tech/badge/bigsmiles)
![license](https://img.shields.io/github/license/dylanwal/bigsmiles)

**(still under development)**

**SMILES** (simplified molecular-input line-entry system) representation is a line notation for molecules with 
given deterministic molecular structures. 

**BigSMILES** is an extension to SMILES which provides support for molecules 
that contain stochastic molecular structures. 

All SMILES are a valid BigSMILES; so this code can also be used with SMILES.

<u>** Features **</u>

* Parse BigSMILES/SMILES into tokens 
* Parse BigSMILES into Compact Graph / Abstract Parse Tree (core data structure for this package)
* Validation of BigSMILES/SMILES (not comprehensive yet)
* Molecular weight distributions
  * Create theory distributions (Poisson, Schulz Zimm, log normal, etc.) 
  * Compute Mn, D, skew, kurtosis, etc.
  * Quick plotting
  * Assign to a stochastic object.
* Generate molecules(SMILES) from a polymer BigSMILES + distribution
* Visualization:
  * Abstract parse tree
  * Compact graph
  
## Documentation

For quickstart, tutorials, reference material, BigSMILES, everything...: 

[**Documention**](https://dylanwal.github.io/BigSMILES/)


---

## Installation

Pip installable package available ([pypi: bigsmiles](https://pypi.org/project/bigsmiles/))

`pip install bigsmiles`


---

## Requirements / Dependencies
Python 3.7 and up

No additional package for core functionality.

Optional dependencies (you will be prompted to install when you use these parts of the code):

* [Numpy](https://github.com/numpy/numpy) (Distributions)
* [Scipy](https://github.com/scipy/scipy) (Distributions)


---

## Basic Usage

The following code generate the core data structure used in this package. During construction the BigSMILES string
will be validated, the Compact graph will be created, and the Abstract parse tree will be generated.

#### Code:
```python
import bigsmiles as bs

polymer_string = "CC{[>][<]CC(C)[>][<]}CC(C)=C"
polymer = bs.BigSMILES(polymer_string)
```


## Features NOT implemented yet
* ladder polymers
* Validation is not comprehensive
* Not all aromatic rings are processed correctly at the moment 
... Lots more to come

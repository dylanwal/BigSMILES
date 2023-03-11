# BigSMILES Parser

---
---

![PyPI](https://img.shields.io/pypi/v/bigsmiles)
![downloads](https://static.pepy.tech/badge/bigsmiles)
![license](https://img.shields.io/github/license/dylanwal/bigsmiles)

**(still under development; but usable)**

**SMILES** (simplified molecular-input line-entry system) representation is a line notation for molecules with 
given deterministic molecular structures. 

**BigSMILES** is an extension to SMILES which provides support for molecules 
that contain stochastic molecular structures. The code here parses the string into and abstract syntax tree.

[**Documention**](https://dylanwal.github.io/BigSMILES/)


---

## Installation

Pip installable package available

`pip install bigsmiles`

[pypi: bigsmiles](https://pypi.org/project/bigsmiles/)


---

## Requirements / Dependencies
Python 3.7 and up


---

## Basic Usage

#### Code:
```python
import bigsmiles as bs

polymer_string = "CC{[>][<]CC(C)[>][<]}CC(C)=C"
polymer = bs.BigSMILES(polymer_string)
```


## Documentation

For quickstart, tutorials, reference material, BigSMILES, everything...  see documentation: 

[**Documention**](https://dylanwal.github.io/BigSMILES/)


## Features NOT implemented yet
* ladder polymers
* Validation is not comprehensive
* Not all aromatic rings are processed correctly at the moment 
... Lots more to come
* 

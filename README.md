# BigSMILES Parser

---
---


**(still under development; but usable)**

SMILES (simplified molecular-input line-entry system) representation is a line notation for molecules with 
given deterministic molecular structures. **BigSMILES** is an extension to SMILES which provides support for molecules 
that contain stochastic molecular structures. The code here parses the string into and abstract syntax tree.

[Documention]()







---

## Installation

Pip installable package available

`pip install bigsmiles`

[pypi: bigsmiles](https://pypi.org/project/bigsmiles/)


---

## Requirements / Dependencies
Python 3.7+


---

## Basic Usage

#### Code:
```python
import bigsmiles as bs

polymer_string = "CC{[>][<]CC(C)[>][<]}CC(C)=C"
polymer = bs.BigSMILES(polymer_string)
```


## Features NOT implemented yet
* ladder polymers
* Validation
  * Is not comprehensive
* Not all aromatic rings are processed correctly at the moment 


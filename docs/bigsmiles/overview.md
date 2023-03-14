---

## Purpose

This specification formally defines the BigSMILES notation (and the specification used by the BigSMILES python package). 

**BigSMILES** is a structurally based line notation based on SMILES (simplified molecular-input line-entry system) 
representation that **adds support for stochasticity**. The primary reason for BigSMILES is to present a
compact yet robust structurally representation that is both machine and human-readable.

BigSMILES is pronounced like: 'big smiles' 

The main motivation for developing BigSMILES was to provide an efficient representation for the  polymer community. 
Polymers are intrinsically stochastic molecules that are often ensembles with a distribution of chemical structures. 
This difficulty limits the applicability of all deterministic representations developed for small molecules. 

This standard was originally published in 
[ACS Cent. Sci. 2019, 5, 9, 1523–1531](https://pubs.acs.org/doi/10.1021/acscentsci.9b00476) and
was updated in this 
[specification](https://olsenlabmit.github.io/BigSMILES/docs/line_notation.html#the-bigsmiles-line-notation).

This documentation builds upon, provides corrections, and further specification to these prior documentations. 




## History / notes

SMILES was originally published in [J. Chem. Inf. Comput. Sci. 1988, 28, 1, 31–36](https://doi.org/10.1021/ci00057a005).
It has since become one of the most popular representations, it however remains imprecise. As a result,
there are dozens of variations and implementations with different interpretations. Among the areas of most ambiguity 
include:

* stereo chemistry (R/S, E/Z)
* aromaticity
* implicit hydrogen counts
* syntax

Among the most prevalent _documented_ revisions include:
* [original publication - Weininger](https://pubs.acs.org/doi/abs/10.1021/ci00057a005)
* [SMILES chapter - Weininger](https://onlinelibrary.wiley.com/doi/epdf/10.1002/9783527618279.ch5)
* [Daylight](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)
* [openSMILES](http://opensmiles.org/opensmiles.html) ([github](https://github.com/opensmiles))
* [SMILES+ - github](https://github.com/IUPAC/IUPAC_SMILES_plus)

Reaction SMILES:
* [Daylight Summer School](https://www.daylight.com/meetings/summerschool01/course/basics/smirks.html)
* [Daylight SMIRKS](https://www.daylight.com/dayhtml/doc/theory/theory.smirks.html)


!!! warning

    Given the variablity in SMILES specifications and implenetations, the results of this python package may
    vary from others. A best effort is made to align with common standards, while supporting the nessesary 
    features of BigSMILES.


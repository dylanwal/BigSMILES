---

## High level

Molecules can be represented at a graph (nodes, edges). Nodes represent atoms, and edges represent bonds. 
In the context of BigSMILES, polymer can similarly be view as a graph.

The construction of a SMILES string involves choosing a starting atom (any atom is a valid starting point) and traverse
the molecule where the traversal path becomes 

![fig1](../img/fig1.svg)


## Atoms

### Privileged Atoms

Privileged atoms is a subset of elements that can be written without square brackets. Members of this group include
alaphatics and aromatic which will be written in lowercase. (aromatics will be discussed more later) 
The atoms have the following properties:

* isotope is taken to be unspecified
* implicit hydrogens will be added to fill the lowest normal valency of the atom
* no stereochemistry 
* atom change is set to zero
* atom class is unspecified

The number of implicit hydrogens is calculated following:

    implicit_hydrogens = valency - number_of_bonds
or

    implicit_hydrogens = valency - sum(bond.bond_order for bond in bonds)


### Extended Atoms

Extended Atoms are written enclosed in square brackets `[]`. 
The Extended Atom pattern consists of 6 parts:

* isotope
* element (required)
* stereo
* hydrogens
* charge
* class

Isotopes is an optional specification place first in the 'Extended Atom' pattern. If not specified, the isotope is 
taken to be the naturally-occurring isotopic ratio. 

Element is a required member of the 'Extended Atom' pattern. It includes the 118 IUPAC defined atomic symbols 
and several aromatic symbols.

Stereochemistry is an optional specification placed after the element symbol. If not present, the stereochemistry
is taken to be undefined. More discussion on stereochemistry is below.

## Bonds

Atoms that are adjacent in a SMILES string are assumed to be joined by a single or aromatic bond in the case of two (see Aromaticity). 
All other bonds must be explicitly written. 

A bond must be proceeded by an Atom.

A single bond can be explicitly represented with '-', but it is rarely necessary.

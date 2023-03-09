
The BigSMILES will be parsed into an abstract syntax tree. There are 3 classes of nodes:

* root node: `BigSMILES`
* intermediate nodes: `StochasticObject`, `StochasticFragment`, `Branch`
* leaf nodes: `BondDescriptorAtom`, `Atom`, `Bond`


BigSMILE objects (`Atom`, `Bond`, `BondingDescriptor`, `Branch`, `StochasticFragment`, `StochasticObject`, `BigSMILES`) 
only holds data. 


!!! info
    Creation of these objects are handled by the 'Constructor' methods.


```mermaid
classDiagram

    class BigSMILES {
        nodes: [intermediate nodes, leaf nodes]
        atoms: [Atom]
        bonds: [Bond]
        rings: [Bond]
    }
    
    
    class StochasticObject {
        nodes: [StochasticFragment]
        bond_left: Bond
        bond_right: Bond
        bond_descriptors: BondDescriptor
    }
    
    
    class StochasticFragment {
        nodes: [intermediate nodes, leaf nodes]
        rings: [Bond]
    }
    
    
    class Branch {
        nodes: [intermediate nodes, leaf nodes]
    }
    
    
    class BondDescriptorAtom {
        BondDescriptor: descriptor
        Bond: bond
    }
    
    
    class BondDescriptor {
        instances: list[BondDescriptorAtom] 
    }

    
    class Bond {
        atom1: Atom
        atom2 Atom
    }
    
    class Atom {
        bonds: list[Bond]: 
    }

    BigSMILES --|> Atom
    BigSMILES --|> Bond
    BigSMILES --|> Branch
    BigSMILES --|> StochasticObject
    StochasticObject --|> StochasticFragment
    StochasticFragment --|> BondDescriptorAtom
    BondDescriptor --|> BondDescriptorAtom
    StochasticFragment --|> Atom
    StochasticFragment --|> Bond
    StochasticFragment --|> Branch
    StochasticFragment --|> StochasticObject
    Branch --|> BondDescriptorAtom
    Branch --|> StochasticObject
    Branch --|> Bond
    Branch --|> Atom
    
```



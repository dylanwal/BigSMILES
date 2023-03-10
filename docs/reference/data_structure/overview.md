
The BigSMILES will be parsed into an abstract syntax tree. There are 3 classes of nodes:

* root node: `BigSMILES`
* intermediate nodes: `StochasticObject`, `StochasticFragment`, `Branch`
* leaf nodes: `BondDescriptorAtom`, `Atom`, `Bond`


BigSMILE objects (`Atom`, `Bond`, `BondingDescriptor`, `Branch`, `StochasticFragment`, `StochasticObject`, `BigSMILES`) 
only holds data. 


!!! info
    Creation of these objects are handled by the 'Constructor' methods.

<b>
The figure shows the data objects and links to other data objects. 
The arrows point from 'parent' objects to 'child' objects. Meaning `Reaction` --> `BigSMILES` means that you will 
find attributes in `Reaction` that are `BigSMILES` objects.

</b>


```mermaid
classDiagram
    
    class Reaction {
        reactants: [BigSMILES]
        agents: [BigSMILES]
        products: [BigSMILES]
    }

    class BigSMILES {
        nodes: [Atom, Bond, Branch, StochasticObject]
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
        nodes: [Atom, Bond, Branch, StochasticObject]
        rings: [Bond]
    }
    
    
    class Branch {
        nodes: [Atom, Bond, Branch, StochasticObject]
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

    Reaction --|> BigSMILES
    BigSMILES --|> Atom
    BigSMILES --|> Bond
    BigSMILES --|> Branch
    BigSMILES --|> StochasticObject
    StochasticObject --|> StochasticFragment
    StochasticObject --|> BondDescriptor
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




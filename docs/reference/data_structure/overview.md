
## Abstract Syntax Tree

root node: `BigSMILES` 

intermediate nodes: `StochasticObject`, `StochasticFragment`, `Branch`

leaf nodes: `BondDescriptorAtom`, `Atom`, `Bond`


BigSMILE objects (`Atom`, `Bond`, `BondingDescriptor`, `Branch`, `StochasticFragment`, `StochasticObject`, `BigSMILES`) 
only holds data. Creation of these objects are handled by the `Constructor`.

```mermaid
classDiagram

    class BigSMILES {
        list: nodes
    }
    
    
    class StochasticObject {
        int: id_
        list: nodes
        BondingDescriptor: end_group_left
        BondingDescriptor: end_group_right
    }
    
    
    class StochasticFragment {
        int: id_
        list: nodes
    }
    
    
    class Branch {
        int: id_
        list: nodes
    }
    
    
    class BondDescriptorAtom {
        int: id_
        BondDescriptor: descriptor
        Bond: bond
    }
    
    
    class BondDescriptor {
        str: symbol
        int: index_
        Enum: type_
        list[BondDescriptorAtom]: instances
    }

    
    class Bond {
        int: id_
        str: symbol
        Enum: type_
        Atom: atom1
        Atom: atom2
        int: ring_id
    }
    
    class Atom {
        int: id_
        str: symbol
        Enum: type_
        int: isotope
        int: charge
        Enum: chiral
        int: valance
        bool: orgainic
        list[Bond]: bonds
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



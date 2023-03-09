### Print Tree
#### Code:
```python
import bigsmiles

polymer_string = "CC{[>][<]CC(C)[>][<]}CC(C)=C"
polymer = bigsmiles.BigSMILES(polymer_string)
polymer.print_tree()
```

#### Output:
```python
BigSMILES: CC{[>][<]CC(C)[>][<]}CC(C)=C
├── Atom: C
├── Bond: 
├── Atom: C
├── Bond: 
├── StochasticObject: {[>][<]CC(C)[>][<]}
│    └── StochasticFragment: [<]CC(C)[>]
│        ├── BondDescriptorAtom: [<]
│        ├── Bond: 
│        ├── Atom: C
│        ├── Bond: 
│        ├── Atom: C
│        ├── Branch: (C)
│        │    ├── Bond: 
│        │    └── Atom: C
│        ├── Bond: 
│        └── BondDescriptorAtom: [>]
├── Bond: 
├── Atom: C
├── Bond: 
├── Atom: C
├── Branch: (C)
│    ├── Bond: 
│    └── Atom: C
├── Bond: =
└── Atom: C
```
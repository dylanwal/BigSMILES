---

::: bigsmiles.methods.display_methods.tree_to_string.tree_to_string


## Example

#### Code:
```python
import bigsmiles
import bigsmiles.methods as bs_methods

polymer_string = "CC{[>][<]CC(C)[>][<]}CC(C)=C"
polymer = bigsmiles.BigSMILES(polymer_string)
text = bs_methods.print_tree(polymer)
print(text)
```

#### Output:
```text
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


## Configuration

::: bigsmiles.methods.display_methods.tree_to_string.TreeConfig
---


!!! info "Tokenization"
    
    Tokenization refers to the process of splitting a input (BigSMILES string) into 
    indivdual terms (Atoms, Bonds, Bonding Descriptors, etc.).


## Text Tokenizer

**Input**: BigSMILES string

**Output**: List of individual tokens as strings. 

```python
import bigsmiles

bigsmiles_str = "CC{[>][<]CC(C)[>][<]}CC(C)=C"

bigsmiles_split_str = bigsmiles.tokenize_text(bigsmiles_str)
print(bigsmiles_split_str)
# output: ['C', 'C', '{', '[>]', '[<]', 'C', 'C', '(', 'C', ')', '[>]', '[<]', '}', 'C', 'C', '(', 'C', ')', '=', 'C']
```

For more information on [bigsmiles.tokenize][bigsmiles.constructors.tokenizer.tokenize_text]

---

## Class Tokenizer

**Input**: BigSMILES string

**Output**: List of Tokens. A [`Token`][[bigsmiles.constructors.tokenizer.Token] has two attributes:

* `Token.kind`: Tokenkind (e.g., `TokenKind.Atom` or `TokenKind.Bond`) 
See [TokenKind][bigsmiles.constructors.tokenizer.TokenKind] (enum) for more information 
* `Token.value`: string (e.g., 'C' or '[>]')

```python
import bigsmiles

bigsmiles_str = "CC{[>][<]CC(C)[>][<]}CC(C)=C"

bigsmiles_tokens = bigsmiles.tokenize(bigsmiles_str)
print(bigsmiles_tokens)
# [
# Atom: C, Atom: C, StochasticStart: {, BondDescriptor: [>], BondDescriptor: [<], Atom: C, Atom: C, 
# BranchStart: (, Atom: C, BranchEnd: ), BondDescriptor: [>], BondDescriptor: [<], StochasticEnd: }, 
# Atom: C, Atom: C, BranchStart: (, Atom: C, BranchEnd: ), Bond: =, Atom: C
# ]

```

For more information on [bigsmiles.tokenize][bigsmiles.constructors.tokenizer.tokenize]

---

## Atom tokenizer

**Input**: BigSMILES string

**Output**: dictionary of attributes 

```python
import bigsmiles

atom_dict = bigsmiles.tokenize_atom_symbol("[13C@H+:1]")
print(atom_dict)
# {'isotope': 13, 'symbol': 'C', 'stereo': '@', 'hydrogens': 1, 'charge': 1, 'class_': 1}
```

For more information on [bigsmiles.tokenize_atom_symbol][bigsmiles.constructors.tokenizer.tokenize_atom_symbol]

---

## Bonding Descriptor Tokenizer

**Input**: Bonding Descriptor string

**Output**: tuple[bond symbol, index]

```python
import bigsmiles
bond_descr = bigsmiles.tokenize_bonding_descriptor("[$1]")
print(bond_descr)
# ('$', 1)
```


For more information on 
[bigsmiles.tokenize_bonding_descriptor][bigsmiles.constructors.tokenizer.tokenize_bonding_descriptor]

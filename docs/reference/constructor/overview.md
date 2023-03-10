---

Constructors are used to build BigSMILES objects. There are three types of the constructor available:

!!! info "Constructor"
    Base constructor that requires full parse information.

    **^^There are two approaches to building a BigSMILES object:^^**

    * build it step-by-step

        1. open branch 
        2. add atoms 
        3. close branch

    * build it in chunks

        1. build the branch in a separate BigSMILES object 
        2. attach the branch at a specified position

    ^^Examples^^

    `add_atom(parent, symbol, isotope, stereo, hydrogens, charge, valence, class_, kwargs)`


!!! info "Constructor String"
    Wraps several of the base constructor constructors. Performs parsing of string inputs before passing values to 
    the base constructor functions.
    
    Methods are appended with '_str' to denote the methods that accept strings.

    ^^Examples^^

    `add_atom(parent, symbol_str, kwargs)`



!!! info "Constructor Tokens"
    Converts a list[Tokens] into a BigSMILES leveraging the 'constructor string'. These function control more of 
    the flow of constructing a BigSMILES; i.e., knowing when to call what constructor function given the situation.

    ^^Example^^
    
    "C(C)C" --> [Token(TokenKind.Atom, "C"), Token(TokenKind.BranchStart, "("), Token(TokenKind.Atom, "C"),
    Token(TokenKind.BranchEnd, ")"), Token(TokenKind.Atom, "C")]

    "[13C@H+:1]" --> {"symbol": "C", "isotope": 13, "stereo": "@", "hydrogens": 1, "charge": 1, "class_": 1}


## Tokenizer

The tokenizer receives BigSMILES strings, Atom strings, or Bonding Descriptor strings and breaks it up into 
the individual tokens or parses in to a dictionary of distinct values. 


## Reaction

The reaction tokenizer parses and tokenized the reaction string. 
Molecules are constructed via `BigSMILES` --> `Constructor Tokens`.

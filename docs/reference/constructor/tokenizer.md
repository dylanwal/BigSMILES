---

## Regular Expression Patterns

Tokenization is mainly done through regular expression patterns.

!!! info "Regular expression Patterns"
    Bond = '-|=|#|$'
    
    Atom = 'S|P|O|N|I|F|Cl|C|Br|B'
    
    Aromatic = 'o|c|n|p|s|b'
    
    AtomExtend = r'(?:\\[)(?P<isotope>[\\d]{1,3})?
    (?P<element>o|c|n|p|s|b|Zr|Zn|Yb|Y|Xe|W|V|U|Ts|Tm|Tl|Ti|Th|Te|Tc|Tb|Ta|Sr
    |Sn|Sm|Si|Sg|Se|Sc|Sb|S|Ru|Rn|Rh|Rg|Rf|Re|Rb|Ra|Pu|Pt|Pr|Po|Pm|Pd|Pb|Pa|P|Os|Og|O|Np|No|Ni|Nh|Ne|Nd|Nb|Na|N|Mt|Mo|Mn|Mg
    |Md|Mc|Lv|Lu|Lr|Li|La|Kr|K|Ir|In|I|Hs|Ho|Hg|Hf|He|H|Ge|Gd|Ga|Fr|Fm|Fl|Fe|F|Eu|Es|Er|Dy|Ds|Db|Cu|Cs|Cr|Co|Cn|Cm|Cl|Cf|Ce
    |Cd|Ca|C|Br|Bk|Bi|Bh|Be|Ba|B|Au|At|As|Ar|Am|Al|Ag|Ac)
    (?P<stereo>@{1,2})?(?P<hydrogens>H[\\d]?)?(?P<charge>[-|\\+]{1,3}[\\d]?)?(?P<class_>:\\d{1,3})?(?:\\])'
    
    BranchStart = r'\('
    
    BranchEnd = r'\)'
    
    Ring = r'[\d]{1}'
    
    Ring2 = r'%[\d]{2}'
    
    BondEZ = r'/|\\'
    
    Disconnected = r"\."
    
    BondDescriptorLadder = r"\[[$<>][\d]\[[$<>][\d]?\][\d]?\]"
    
    BondDescriptor = r"\[[$<>][\d]?[\d]?\]"
    
    StochasticSeperator = r",|;"
    
    StochasticStart = r'\{'
    
    StochasticEnd = r'\}'
    
    ImplictEndGroup = r'\[\]'
    
    Rxn = r'>>|>'




## Atom tokenizer

::: bigsmiles.constructors.tokenizer.tokenize_atom_symbol


---

## Bonding Descriptor Tokenizer

::: bigsmiles.constructors.tokenizer.tokenize_bonding_descriptor


---

## BigSMILES Tokenizer

The tokenizer leverage python's regular expression to identify and label valid BigSMILE tokens.


::: bigsmiles.constructors.tokenizer.tokenize


::: bigsmiles.constructors.tokenizer.Token


::: bigsmiles.constructors.tokenizer.TokenKind
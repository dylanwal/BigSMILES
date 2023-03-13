

## Grammar Notation

SMILES can be contextualized as a [context-free grammar](https://en.wikipedia.org/wiki/Context-free_grammar).
This means that the syntax can be  defined by:

`A --> <span>&#124;</span>`

[extended Backus–Naur form (EBNF)](https://en.wikipedia.org/wiki/Extended_Backus%E2%80%93Naur_form) is a notation 
for formally describing syntax: how to write the linguistic We will use EBNF to describe the
features in a language. Using this notation, a person can determine
whether a program is syntactically correct: whether it adheres to the grammar
and punctuation rules of the language.

An EBNF description is an unordered list of EBNF rules. Each EBNF rule has three parts: 

 * a left–hand side (LHS); names the EBNF rule
 * a right-hand side (RHS); supplies a description
 * the ::= character separating these two sides; read this symbol as “is defined as”.



The following will provide a brief definition of the syntax used.

A valid 'word' in a grammar is defined as containing only terminal symbols. In the formulation of EBNF a terminal is 
a string and will only appear in the right-hand side of a rule.
but never appears on the left-hand side in the main grammar, 
although it may appear on the left-hand side of a rule in the grammar for terminals.


` symbol ::= expression`


| Syntax   | interpretation                                                          |
|----------|-------------------------------------------------------------------------|
| "string" | matches the sequence of characters that appear inside the double quotes |
| 'string' | matches the sequence of characters that appear inside the single quotes |
| A B   | matches A followed by B                                                 |
| A  <span>&#124;</span> B | matches A or B but not both                                             |
| A+ | matches one or more occurrences of A                                    |
| A* | matches zero or more occurrences of A                                   |
| A? | matches A or nothing; optional A                                        |
|(A) | A is treated as a unit and may be combined as described in this list    |


To prove that a symbol is legal according to some EBNF rule, we must match
all its characters with all the items in the EBNF rule, according to that rule’s
description. If there is an exact match —we exhaust the characters in the
symbol at the same time when exhaust the rule’s description— we classify
the symbol as legal according to that EBNF description and say it matches;
otherwise we classify the symbol as illegal and say it doesn’t match.




## General



```text
terminator ::= SPACE | TAB | LINEFEED | CARRIAGE_RETURN | END_OF_STRING


digit_not_zero_one ::= '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9'
digit_not_zero ::= '1' | digit_not_zero_one
digit ::= '0' | digit_not_zero
two_digit ::= digit_not_zero digit
three_digit ::= digit_not_zero digit digit

aliphatic ::= 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I'
aromatic ::= 'b' | 'c' | 'n' | 'o' | 's' | 'p'
extended_atom ::= '[' isotope? symbol chiral? hydrogens? charge? class? ']'
isotope ::= digit_not_zero | two_digit | three_digit
element_symbols ::= 'H' | 'He' | ... | 'Ts' | 'Og'
chiral ::= '@' | '@@'
hydrogens ::= 'H' | 'H' digit_not_zero_one
charge ::= '-' | '-' digit_not_zero_one |'+' | '+' digit_not_zero_one |'--' | ... | '---------' | '++' | ... | '+++++++++'
class ::= ':' digit | ':' two_digit | ':' three_digit

atom ::= aliphatic | aromatic | extended_atom | '*'


bond ::= '-' | '=' | '#' | '$' | ':' | '/' | '\'
disconnect ::= '.'


branch ::= '(' chain ')' | '(' bond chain ')' | '(' disconnect chain ')'
ring ::= bond? digit_not_zero | bond? '%' two_digit | bond? '%' '('three_digit ')'

branched_atom ::= atom ring* branch*


bonding_descriptor_symbol ::= '<' | '>' | '$'
bonging_descriptor_symbol_index ::= '[' bonding_descriptor_symbol (digit_not_zero | two_digit)? ']'


stochastic_fragment ::= chain? (bond | disconnect)? bonding_descriptor chain? (bond | disconnect)? bonding_descriptor?
stochastic_object ::= '{' bonding_descriptor stochastic_fragment (',' stochastic_fragment)? bonding_descriptor '}'

chain ::= branched_atom | chain branched_atom | chain bond branched_atom | chain dot branched_atom | stochastic_object | chain stochastic_object | chain bond stochastic_object
stochastic_object_implicit ::= '{[]' stochastic_fragment (',' stochastic_fragment)? '[]}'
bigSMILES ::= chain+ terminator | stochastic_object_implicit terminator


reaction ::= bigSMILES (',' BigSMILES)? '>>' bigSMILES (',' BigSMILES)? | bigSMILES (',' BigSMILES)? '>' bigSMILES (',' BigSMILES)? '>' bigSMILES (',' BigSMILES)?
```


tabular proof and derivation trees
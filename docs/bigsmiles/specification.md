---



## Atoms

### Privileged Atoms

<table id="privileged_atoms">
    <tr>
        <th> class </th>
        <th> grammar </th>
    </tr>
    <tr>
        <th> aliphatic </th>
        <th>'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I'</th>
    </tr>
    <tr>
        <th> aromatic </th>
        <th>'b' | 'c' | 'n' | 'o' | 's' | 'p'</th>
    </tr>
</table>

### Extended Atoms

`'[' isotope? element chiral? hydrogens? charge? class? ']'`

<table id="extended_atoms">
    <tr>
        <th> class </th>
        <th> grammar </th>
    </tr>
    <tr>
        <th> isotope </th>
        <th> '[0-9]{0,3}' </th>
    </tr>
    <tr>
        <th> element </th>
        <th> 
"o' | 'n' | 'c' | 'p' | 'b' | 's' | 'Zr' | 'Zn' | 'Yb' | 'Y' | 'Xe' | 'W' | 'V' | 'U' | 'Ts' | 'Tm' | 'Tl' | 'Ti' | 
'Th' | 'Te' | 'Tc' | 'Tb' | 'Ta' | 'Sr' | 'Sn' | 'Sm' | 'Si' | 'Sg' | 'Se' | 'Sc' | 'Sb' | 'S' | 'Ru' | 'Rn' | 'Rh' | 
'Rg' | 'Rf' | 'Re' | 'Rb' | 'Ra' | 'Pu' | 'Pt' | 'Pr' | 'Po' | 'Pm' | 'Pd' | 'Pb' | 'Pa' | 'P' | 'Os' | 'Og' | 'O' | 
'Np' | 'No' | 'Ni' | 'Nh' | 'Ne' | 'Nd' | 'Nb' | 'Na' | 'N' | 'Mt' | 'Mo' | 'Mn' | 'Mg' | 'Md' | 'Mc' | 'Lv' | 'Lu' | 
'Lr' | 'Li' | 'La' | 'Kr' | 'K' | 'Ir' | 'In' | 'I' | 'Hs' | 'Ho' | 'Hg' | 'Hf' | 'He' | 'H' | 'Ge' | 'Gd' | 'Ga' | 
'Fr' | 'Fm' | 'Fl' | 'Fe' | 'F' | 'Eu' | 'Es' | 'Er' | 'Dy' | 'Ds' | 'Db' | 'Cu' | 'Cs' | 'Cr' | 'Co' | 'Cn' | 'Cm' | 
'Cl' | 'Cf' | 'Ce' | 'Cd' | 'Ca' | 'C' | 'Br' | 'Bk' | 'Bi' | 'Bh' | 'Be' | 'Ba' | 'B' | 'Au' | 'At' | 'As' | 'Ar' | 
'Am' | 'Al' | 'Ag' | 'Ac" 
</th>
    </tr>
    <tr>
        <th> chiral </th>
        <th> '@' | '@@' </th>
    </tr>
    <tr>
        <th> hydrogens </th>
        <th> 'H{0,3}' </th>
    </tr>
    <tr>
        <th> charge </th>
        <th> '[-|+]{1,3}' | '[-|+][0-9]' </th>
    </tr>
    <tr>
        <th> class </th>
        <th> ':[0-9]{1,3}' </th>
    </tr>

</table>


## Bonds

<table id="bond">
<tr>
        <th colspan="2"> <b>Bonds</b> </th>
    </tr>
    <tr>
        <th> Bond </th>
        <th>'-' | '=' | '#' | '$' | ':' | '/' | '\' </th>
    </tr>
    <tr>
        <th> Ring </th>
        <th>'[1-9]' |  '%[1-9][0-9]' </th>
    </tr>
    <tr>
        <th> Disconnection </th>
        <th> '.' </th>
    </tr>
</table>

### Branch

<table>
    <tr>
        <th colspan="2"> <b>Branch</b> </th>
    </tr>
    <tr>
        <th colspan="2"> <b> Compound Elements </b> </th>
    </tr>
    <tr>
        <th> BigSMILES </th>
        <th> aliphatic | aromatic | Extended Atom | Bond</th>
    </tr>
    <tr>
        <th colspan="2"> <b> BigSMILES </b> </th>
    </tr>
    <tr>
        <th> BigSMILES </th>
        <th> aliphatic | aromatic | Extended Atom | Bond</th>
    </tr>
    <tr>
        <th> BigSMILES </th>
        <th> terminator | chain terminator </th>
    </tr>
    <tr>
        <th colspan="2"> <b> BigSMILES Reaction </b> </th>
    </tr>
    <tr>
        <th> agent free reaction </th>
        <th> BigSMILE (',' '>>' </th>
    </tr>
</table>






<style> 
 
    #atom {
      border: 2px solid;
    }
    
    #atom tr:hover {background-color: #ddd;}
    
    #atom th {
      padding-top: 12px;
      padding-bottom: 12px;
      text-align: left;
    }

</style>
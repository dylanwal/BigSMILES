import bigsmiles

reaction_string = "C=Cc1ccccc1.C[CH-](.[Li+])CC>Cc1ccccc1>CC(CC){[>][<]CC(c1ccccc1)[>][<]}[H]"
rxn = bigsmiles.Reaction(reaction_string)

print([str(chem) for chem in rxn.reactants])
print([str(chem) for chem in rxn.agents])
print([str(chem) for chem in rxn.products])

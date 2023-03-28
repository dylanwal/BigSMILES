from __future__ import annotations
import copy

import numpy as np

from bigsmiles.data_structures.bigsmiles import BigSMILES, StochasticObject
import bigsmiles.constructors.constructor as constructor
from bigsmiles.data_structures.polymer import Polymer, StochasticObjectSpecification


class GenerationData:
    """ grouping of data for generating a molecules from a Polymer. """
    def __init__(self, polymer: Polymer, rng: np.random.Generator = None):
        self.polymer = polymer
        self.rng = rng
        self.sorted_spec = sorted(polymer.spec, key=lambda x: x.depth)
        self.stoch_depth_level_one_index = self.get_stoch_objects_index()

        # single generation data (re-written each generation)
        self.replacement_fragments: list[BigSMILES] = []

    def get_stoch_objects_index(self):
        return [i for i, obj in enumerate(self.polymer.bigsmiles.nodes) if isinstance(obj, StochasticObject)]

    def generate_stochastic_objects(self):
        """ loop through each stochastic object and generate specific instances of each """
        # use sorted list so innermost stochastic objects are built first and can be used by following
        # parent stochastic object
        for spec_ in self.sorted_spec:
            self.build_fragment(spec_)

    def build_fragment(self, spec: StochasticObjectSpecification):
        """ build a specific realization of a stochastic object. """
        # TODO: limit homo-linear polymer
        N = spec._get_N_for_generation(self.rng)

        bigsmiles_ = BigSMILES()
        for i in range(N):
            chunk = self.get_chunk(spec, i)
            constructor.append_bigsmiles_fragment(bigsmiles_, chunk, bond_symbol="")

        self.replacement_fragments.append(bigsmiles_)

    def get_chunk(self, spec, i) -> BigSMILES:
        if len(spec.stochastic_fragments) == 1:
            return self.build_chunk(stochastic_fragment, i)


def generate_molecules(polymer: Polymer, n: int = 1, rng: np.random.Generator = None) -> BigSMILES | list[BigSMILES]:
    """
    Generate molecules. In this context, a molecule is a specific realization of a polymer. It is calculated by using
    molecular weight distribution and composition information.

    Parameters
    ----------
    polymer:
        polymer you want to generate specific realization for
    n:
        how many molecules you want generated
    rng:
        random number generator

    Returns
    -------
    bigsmiles:
        generated molecules

    """
    data = GenerationData(polymer, rng)

    if n == 1:
        return generate_one_molecule(data)

    molecules = []
    for i in range(n):
        molecules.append(generate_one_molecule(data))

    return molecules


def generate_one_molecule(data: GenerationData) -> BigSMILES:
    """ Generate a single molecule. """
    # build stochastic objects
    data.generate_stochastic_objects()

    # build molecule
    bigsmiles = copy.deepcopy(data.polymer.bigsmiles())

    # loop through top level and replace all stochastic objects with specific realization
    for index, fragment in zip(data.stoch_depth_level_one_index, data.replacement_fragments):
        constructor.replace_stochastic_object(bigsmiles.nodes[index], fragment)

    return bigsmiles

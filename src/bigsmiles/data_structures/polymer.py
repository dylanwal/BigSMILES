from __future__ import annotations
import statistics


from bigsmiles.data_structures.bigsmiles import BigSMILES, StochasticObject, StochasticFragment
from bigsmiles.distributions.base import Distribution


def extract_text_between_pipe(string: str) -> tuple[list[str, ...], list[str, ...]]:
    removed_text = []
    remaining_text = []
    while "|" in string:
        start = string.find("|")
        end = string.find("|", start + 1)
        removed_text.append(string[start + 1:end])
        remaining_text.append(string[:start] + string[end + 1:])
        string = string[:start] + string[end + 1:]
    remaining_text.append(string)

    return removed_text, remaining_text


class StochasticObjectSpecification:
    """ Object contains molecular weight and composition data for a single stochastic object. """
    DEFAULT_N = 10

    def __init__(self, stochastic_object: StochasticObject):
        """
        Parameters
        ----------
        stochastic_object:
            stochastic object
        """
        self.stochastic_object = stochastic_object
        self.depth = self.stochastic_object.nesting_depth_stochastic_objects

        self._N: int | None = None
        self._distribution: int | None = None
        self._composition: float | None = None
        self._repeat_MW: float | None = None

    def __str__(self):
        return str(self.stochastic_object)

    @property
    def stochastic_fragments(self) -> list[StochasticFragment]:
        return self.stochastic_object.nodes

    @property
    def N(self) -> int | None:
        """ If None, it will be obtained from distribution, if distrtion is None, N is default value. """
        return self._N

    @N.setter
    def N(self, N: int):
        if isinstance(N, int) and N > 0:
            self._N = N
        else:
            raise ValueError("Invalid N.")

    @property
    def distribution(self) -> Distribution | None:
        """ molecular weight distribution """
        return self._distribution

    @distribution.setter
    def distribution(self, distribution: Distribution):
        if isinstance(distribution, Distribution):
            self._distribution = distribution
            if self._distribution.repeat_MW is None:
                self._distribution.repeat_MW = self.get_repeat_MW()
        else:
            raise ValueError("Invalid distribution.")

    @property
    def composition(self) -> float | None:
        return self._composition

    @composition.setter
    def composition(self, composition: float):
        self._composition = composition

    def get_repeat_MW(self):
        if self.composition is None:
            return statistics.fmean([fragment.molar_mass for fragment in self.stochastic_object.nodes])

    def _get_N_for_generation(self, rng=None) -> int:
        """

        Parameters
        ----------
        rng: np.random.Generator

        Returns
        -------

        """
        if self.N is not None:
            return self.N
        if self.distribution is not None:
            return self.distribution.draw_N(rng)

        return self.DEFAULT_N


def get_all_stochastic_objects(obj, stochastic_objects: list[StochasticObject] = None) \
        -> list[StochasticObjectSpecification]:
    if stochastic_objects is None:
        stochastic_objects = []

    for node in obj.nodes:
        if isinstance(node, StochasticObject):
            stochastic_objects.append(StochasticObjectSpecification(node))
        if hasattr(node, 'nodes'):
            stochastic_objects += get_all_stochastic_objects(node, stochastic_objects)  # recursive

    return stochastic_objects


class Polymer:
    """ BigSMILES + molecular weight and composition information """
    def __init__(self, text: str = None):
        """
        Parameters
        ----------
        text:
            BigSMILES string
        """
        if '|' in text:
            bigsmiles_text, gen_text = extract_text_between_pipe(text)
            text = "".join(bigsmiles_text)

        self.bigsmiles = BigSMILES(text)
        self.spec: list[StochasticObjectSpecification] = get_all_stochastic_objects(self.bigsmiles)

        # add BigSMILES extensions into spec

    def __str__(self):
        return str(self.bigsmiles)

    def __repr__(self):
        return str(self.bigsmiles)

    # @property
    # def stochastic_objects(self) -> list[StochasticObject]:
    #     return [stoch_obj.stochastic_object for stoch_obj in self._stochastic_objects_extend]

    def get_BigSMILES_gen_str(self) -> str:
        pass

import gillespy2 as gillespy
import numpy as np
import math

import base

# test with number of species = 5, 10, 20


class BIRTH_DEATH(base.BaseCRNModel):
    """Class for BIRTH_DEATH network from the DeepCME paper."""

    params = {
        'kappa': 10.,
        'gamma': 1.
    }

    def __init__(self, endtime, timestep):
        """
        Initialize the model.

        Parameters
        ----------
        endtime : endtime of simulations
        timestep : time-step of simulations
        """
        super().__init__(
            endtime=endtime,
            timestep=timestep,
            model_name="SIR",
        )

        kappa = gillespy.Parameter(name='kappa', expression=self.params['kappa'])
        gamma = gillespy.Parameter(name='gamma', expression=self.params['gamma'])
        self.add_parameter([kappa, gamma])

        # Species
        S1_start = gillespy.Species(name='S1_start', initial_value=1)
        S1_alive = gillespy.Species(name='S1_alive', initial_value=1)
        S1_dead = gillespy.Species(name='S1_dead', initial_value=1)

        S2_start = gillespy.Species(name='S2_start', initial_value=1)
        S2_alive = gillespy.Species(name='S2_alive', initial_value=1)
        S2_dead = gillespy.Species(name='S2_dead', initial_value=1)

        S3_start = gillespy.Species(name='S3_start', initial_value=1)
        S3_alive = gillespy.Species(name='S3_alive', initial_value=1)
        S3_dead = gillespy.Species(name='S3_dead', initial_value=1)

        S4_start = gillespy.Species(name='S4_start', initial_value=1)
        S4_alive = gillespy.Species(name='S4_alive', initial_value=1)
        S4_dead = gillespy.Species(name='S4_dead', initial_value=1)

        S5_start = gillespy.Species(name='S5_start', initial_value=1)
        S5_alive = gillespy.Species(name='S5_alive', initial_value=1)
        S5_dead = gillespy.Species(name='S5_dead', initial_value=1)
        self.add_species([
                S1_start, S1_alive, S1_dead,
                S2_start, S2_alive, S2_dead,
                S3_start, S3_alive, S3_dead,
                S4_start, S4_alive, S4_dead,
                S5_start, S5_alive, S5_dead
            ])

        # Reactions
        # births
        birth_s1 = gillespy.Reaction(
            name='birth_s1',
            reactants={S1_start: 1},
            products={S1_alive: 1},
            rate=kappa
        )
        birth_s2 = gillespy.Reaction(
            name='birth_s2',
            reactants={S2_start: 1},
            products={S2_alive: 1},
            rate=kappa
        )
        birth_s3 = gillespy.Reaction(
            name='birth_s3',
            reactants={S3_start: 1},
            products={S3_alive: 1},
            rate=kappa
        )
        birth_s4 = gillespy.Reaction(
            name='birth_s4',
            reactants={S4_start: 1},
            products={S4_alive: 1},
            rate=kappa
        )
        birth_s5 = gillespy.Reaction(
            name='birth_s5',
            reactants={S5_start: 1},
            products={S5_alive: 1},
            rate=kappa
        )

        # deaths
        death_s1 = gillespy.Reaction(
            name='death_s1',
            reactants={S1_alive: 1},
            products={S1_dead: 1},
            rate=gamma
        )
        death_s2 = gillespy.Reaction(
            name='death_s2',
            reactants={S2_alive: 1},
            products={S2_dead: 1},
            rate=gamma
        )
        death_s3 = gillespy.Reaction(
            name='death_s3',
            reactants={S3_alive: 1},
            products={S3_dead: 1},
            rate=gamma
        )
        death_s4 = gillespy.Reaction(
            name='death_s4',
            reactants={S4_alive: 1},
            products={S4_dead: 1},
            rate=gamma
        )
        death_s5 = gillespy.Reaction(
            name='death_s5',
            reactants={S5_alive: 1},
            products={S5_dead: 1},
            rate=gamma
        )

        
        self.add_reaction([
            birth_s1, birth_s2, birth_s3, birth_s4, birth_s5,
            death_s1, death_s2, death_s3, death_s4, death_s5
        ])
        nb_of_steps = int(math.ceil((endtime / timestep))) + 1
        self.timespan(np.linspace(0, endtime, nb_of_steps))

    def set_species_initial_value(self, species_initial_value):
        """
        Set initial values to species.

        Parameters
        ----------
        species_initial_value : list or 1d array of values, size should be equal
            to the number of species, the order should be coherent with theo order
            of species returned by get_initial_state method.

        Returns
        -------
        None

        """
        self.listOfSpecies['S1_start'].initial_value = species_initial_value[0]
        self.listOfSpecies['S1_alive'].initial_value = species_initial_value[1]
        self.listOfSpecies['S1_dead'].initial_value = species_initial_value[2]

        self.listOfSpecies['S2_start'].initial_value = species_initial_value[3]
        self.listOfSpecies['S2_alive'].initial_value = species_initial_value[4]
        self.listOfSpecies['S2_dead'].initial_value = species_initial_value[5]

        self.listOfSpecies['S3_start'].initial_value = species_initial_value[6]
        self.listOfSpecies['S3_alive'].initial_value = species_initial_value[7]
        self.listOfSpecies['S3_dead'].initial_value = species_initial_value[8]

        self.listOfSpecies['S4_start'].initial_value = species_initial_value[9]
        self.listOfSpecies['S4_alive'].initial_value = species_initial_value[10]
        self.listOfSpecies['S4_dead'].initial_value = species_initial_value[11]

        self.listOfSpecies['S5_start'].initial_value = species_initial_value[12]
        self.listOfSpecies['S5_alive'].initial_value = species_initial_value[13]
        self.listOfSpecies['S5_dead'].initial_value = species_initial_value[14]

    @staticmethod
    def get_species_names():
        """
        Returns list of species names.

        Returns
        -------
        list of all species names. The order of names should be coherent
        with the list returned by get_initial_state method.

        """
        return [
            'S1_start', 'S1_alive', 'S1_dead',
            'S2_start', 'S2_alive', 'S2_dead',
            'S3_start', 'S3_alive', 'S3_dead',
            'S4_start', 'S4_alive', 'S4_dead',
            'S5_start', 'S5_alive', 'S5_dead'
        ]

    @staticmethod
    def get_initial_state():
        pass

    @classmethod
    def get_n_species(cls):
        """Total number of species."""
        return len(cls.get_species_names())

    @classmethod
    def get_initial_settings(cls, n_settings, sigm=None):
        """
        Generate a set of (random) initial states.

        Parameters
        ----------
        n_settings : number of settings to produce.
        sigm : not used

        Returns
        -------
        array of n_settings initial states

        """
        n_species = cls.get_n_species()
        settings = None
        for i in range(n_settings):
            curr_settings = []
            for j in range(n_species):
                if j % 3 == 0:
                    curr_settings.append(np.random.randint(low=10, high=100))
                else:
                    curr_settings.append(0)
            curr_settings = np.array(curr_settings)
            if i == 0:
                settings = curr_settings
            else:
                settings = np.vstack([settings, curr_settings])

        # settings = np.random.randint(low=30, high=200, size=(n_settings, n_species))
        # settings = np.zeros()
        settings = np.reshape(settings, (n_settings, n_species))
        return settings

    @classmethod
    def get_histogram_bounds(cls, species_names_list=None):
        """
        Returns bounds for species histograms.

        Parameters
        ----------
        species_names_list: not used

        Returns
        -------
        histogram_bounds: list of [min, max] values for species

        """
        n_species_for_histogram = len(cls.get_species_for_histogram())
        histogram_bounds = [[0.5, 200.5]] * n_species_for_histogram
        return histogram_bounds

    @staticmethod
    def get_species_for_histogram():
        """
        Returns subset of species of interest.

        Returns
        -------
        list of species names.

        """
        return [
            'S1_alive', 'S1_dead',
            'S2_alive', 'S2_dead',
            'S3_alive', 'S3_dead',
            'S4_alive', 'S4_dead',
            'S5_alive', 'S5_dead'
        ]

    # @classmethod
    # def get_randomized_parameters(cls, param_names, n_settings, sigm=0.5):
    #     randomized = {}
    #     for name in param_names:
    #         if name not in cls.params:
    #             raise KeyError(f"Could not find param {name} in {cls.__name__} class `params` dict.")
    #         val = float(cls.params[name])
    #
    #         randomized[name] = np.random.uniform(val * (1. - sigm), val * (1. + sigm), n_settings)
    #     return randomized
    #
    # def set_parameters(self, params_dict):
    #     for name, val in params_dict.items():
    #         if name not in self.listOfParameters:
    #             raise KeyError(
    #                 f"Could not find {name} parameter in {self.__class__.__name__} model listOfParameters.")
    #         self.set_parameter(name, str(val))

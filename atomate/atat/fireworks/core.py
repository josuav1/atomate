# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

try:
    from math import gcd
except ImportError:
    from fractions import gcd

import multiprocessing
import os
from fireworks import Firework
from fireworks import ScriptTask
from pymatgen.analysis.structure_matcher import StructureMatcher, OccupancyComparator
from monty.serialization import dumpfn
from pymatgen.transformations.standard_transformations import DiscretizeOccupanciesTransformation
from fractions import Fraction

__author__ = 'Matthew Horton'
__credits__ = 'Josua Vieten'
__email__ = 'mkhorton@lbl.gov'

class McsqsFW(Firework):

    def __init__(self, disordered_struct,
                 name="mcsqs",
                 size=None,
                 clusters=None,
                 user_input_settings=None,
                 max_denominator=8,
                 walltime=2.5,
                 ncores=None,
                 **kwargs):
        """
        Find a best SQS approximation for a given disordered
        structure using the mcsqs code from the ATAT library.

        One unusual feature of this Firework: walltime must
        be specified. This is because there are often cases
        where mcsqs will never terminate (it will keep trying
        structures at random); however, whatever structure
        it finds in this time (the 'bestsqs') is often
        still useful.

        :param disordered_struct: disordered structure
        :param name: name of the Firework
        :param size: number of atoms in the output supercell,
        if None will choose the smallest possible supercell
        :param clusters: dict of cluster sizes, e.g.
        {2: 1.5, 3: 1.0} to search for pairs of atoms
        within a radius of 1.5 Å and triplets of atoms
        within a radius of 1.0 Å, if None will try to
        pick sensible defaults, the value of mcsqs is
        that this should converge rapidly with number of
        clusters (aim for around a dozen)
        :param user_input_settings: dict of additional keyword
        value pairs to pass to mcsqs, such as Monte Carlo
        temperature, no validation is performed on this dict
        :param: max_denominator (int):
        Will round occupancies to the nearest rational number,
        with the maximum denominator as specified. Use None not
        to do perform rounding, though not mcsqs might not like this!
        :param: walltime (int): Time to run mcsqs for in minutes,
        2 minutes will be deducted from this to run the database
        tasks.
        :param: ncores (int): number of instances of mcsqs
        to run (mcsqs is not parallelized, these are run
        independently of each other), by default will try
        to detect automatically
        :param kwargs: Other kwargs that are passed to Firework.__init__.
        """

        if walltime < 3:
            raise ValueError("Please increase walltime.")

        if disordered_struct.is_ordered:
            raise ValueError("You must input a disordered structure.")

        if max_denominator:
            trans = DiscretizeOccupanciesTransformation(max_denominator)
            disordered_struct = trans.apply_transformation(disordered_struct)

        if size is None:
            size = self._determine_min_cell(disordered_struct)

        # rescale lattice for unit volume
        disordered_struct.scale_lattice(1.0)

        if clusters is None:
            # TODO: make cluster determination a bit smarter!
            lattice_param = min(disordered_struct.lattice.abc)
            clusters = {
                -2: lattice_param*1.501,
                -3: lattice_param*1.501,
                -4: lattice_param*1.001
            }
        else:
            for cluster in clusters.keys():
                if cluster not in [2, 3, 4, 5, 6]:
                    raise ValueError("Mcsqs only supports clusters of 2-6 atoms.")

        disordered_struct.to(filename='rndstr.in')

        cluster_str = " ".join(["{}={}".format(k, v) for k, v in clusters.items()])
        generate_cluster_cmd = "mcsqs {}".format(cluster_str)

        if user_input_settings:
            user_input_settings_str = " ".join(["{}={}".format(k, v) for k, v
                                            in user_input_settings.items()])
        else:
            user_input_settings_str = ""

        # using same approach as Custodian to detect ncores
        ncores = os.environ.get('NSLOTS') or multiprocessing.cpu_count()
        ncores = int(ncores)

        # mcsqs is not parallelised, and is not expected to finish within the walltime
        # it's a Monte Carlo process that can run indefinitely, so we use timeout
        # to ensure the command finishes

        run_mcsqs_cmd = """
for (( id=0 ; id<{} ; id ++ ))

do
    timeout {}m mcsqs -n {} {}  -ip=$id & 
done
""".format(ncores, walltime - 2, size, user_input_settings_str)

        # command to find the best SQS among the several instances of mcsqs run in parallel
        get_bestsqs_cmd = "mcsqs -best"

        # write the mcsqs version to a file for provenance
        write_version_cmd = "mcsqs -v | head -n 1 > mcsqs_version.txt"

        # write our input args, so that they can be picked up by the drone
        dumpfn({
            'clusters': clusters,
            'user_input_settings': user_input_settings,
            'walltime': ncores*(walltime - 2)
        }, "mcsqs_input_args.json")

        tasks = [
            ScriptTask(script=[
                run_mcsqs_cmd,
                get_bestsqs_cmd,
                write_version_cmd,
            ], shell_exe='/bin/bash'),
        ]

        super(McsqsFW, self).__init__(tasks, name=name, **kwargs)

    @staticmethod
    def _determine_min_cell(disordered_struct):
        """
        Get a minimum cell size for an ordered structure.
        """
        denoms = {}
        for sp_occu in disordered_struct.species_and_occu:
            for sp, occu in sp_occu.items():
                denom = Fraction(occu).limit_denominator(100).denominator
                if sp in denoms:
                    denoms[sp] = max(denom, denoms[sp])
                else:
                    denoms[sp] = denom
        return max(denoms.values())
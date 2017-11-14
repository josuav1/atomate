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
from atomate.atat.firetasks.parse_outputs import McsqsToDbTask

__author__ = 'Matthew Horton'
__credits__ = 'Josua Vieten'
__email__ = 'mkhorton@lbl.gov'

class DbtaskFW(Firework):

    def __init__(self,path,name="mcsqs", 
                 **kwargs):
        """
        Save results to database.

        :param kwargs: Other kwargs that are passed to Firework.__init__.

        """

        db_fw = McsqsToDbTask(calc_dir=path,db_file='sqsdb.json')

        tasks = [db_fw]

        super(DbtaskFW, self).__init__(tasks, name=name, **kwargs)

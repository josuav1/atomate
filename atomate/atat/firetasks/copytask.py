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

class CopytaskFW(Firework):

    def __init__(self,target,name="mcsqs", 
                 **kwargs):
        """
        Copy resulting bestsqs.out, clusters.out, rndstr.in, bestcorr.out, rndstrgrp.out, and sqscell.out to working sub-directory

        target: target directory

        :param kwargs: Other kwargs that are passed to Firework.__init__.

        """

        filelist = ['bestsqs.out', 'clusters.out', 'rndstr.in', 'bestcorr.out', 'rndstrgrp.out', 'sqscell.out', 'mcsqs_version.txt', 'mcsqs_input_args.json']

        taskstring = ''
        for i in range (0,len(filelist)):
            taskstring = taskstring + "cp -f ./" + filelist[i] + " " + target + " ;"

        tasks = [ScriptTask(script=[taskstring], shell_exe='/bin/bash')]

        super(CopytaskFW, self).__init__(tasks, name=name, **kwargs)

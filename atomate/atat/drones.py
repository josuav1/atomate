# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import datetime
import os

from pymatgen.core.structure import Structure
from atomate.utils.utils import get_logger
from atomate import __version__ as atomate_version
from monty.serialization import loadfn

__author__ = 'Matthew Horton'
__credits__ = 'Kiran Mathew, Shyue Ping Ong, Shyam Dwaraknath, Anubhav Jain, Josua Vieten'
__email__ = 'mkhorton@lbl.gov'
# based on VaspDrone

logger = get_logger(__name__)


class McsqsDrone:

    def assimilate(self, path):

        logger.info("Assimilating mcsqs for base dir: {}".format(path))

        input_structure = Structure.from_file(os.path.join(path, 'rndstr.in'))
        anonymous_formula = input_structure.composition.anonymized_formula

        output_structure = Structure.from_file(os.path.join(path, 'bestsqs.out'))

        # get our objective function
        with open(os.path.join(path, 'bestcorr.out')) as f:
            lines = f.read().split('\n')
            lines = [l for l in lines if l]
            try:
                obj_func = float(lines[-1].split('=')[1])
            except:
                obj_func = str(lines[-1].split('=')[1])

        # get scaling matrix
        try:
            l = input_structure.lattice
            l_sqs = output_structure.lattice
            scaling_matrix = l.find_mapping(l_sqs)[2].tolist()
        except:
            scaling_matrix = "Could not determine scaling matrix?"

        # get mcsqs version
        try:
            # if running via atomate, we will create this file
            with open(str(path) +'mcsqs_version.txt') as f:
                mcsqs_version = str(f.read())
        except:
            mcsqs_version = "Unknown mcsqs version"

        # get number of clusters
        with open(str(path) +'clusters.out') as f:
            lines = f.read().split('\n')
            num_clusters = (len(lines) - 2) / 7

        # load input args
        try:
            input_args = loadfn(os.path.join(path, "mcsqs_input_args.json"))
            walltime = input_args['walltime']
            clusters = input_args['clusters']
            user_input_settings = input_args['user_input_settings']
        except:
            walltime = None
            clusters = None
            user_input_settings = None

        return {
            'anonymous_formula': anonymous_formula,
            'disordered': input_structure.as_dict(),
            'bestsqs': output_structure,
            'clusters': clusters,
            'num_clusters': num_clusters,
            'user_input_settings': user_input_settings,
            'objective_function': obj_func,
            'walltime': walltime,
            'atomate_version': atomate_version,
            'mcsqs_version': mcsqs_version,
            'spacegroup': input_structure.get_space_group_info()[0],
            'scaling_matrix': scaling_matrix,
            'size': len(output_structure)/len(input_structure),
            'last_updated': str(datetime.datetime.now())
        }

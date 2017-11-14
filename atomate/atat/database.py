# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines the database classes.
"""
import datetime
import os
from atomate.utils.database import CalcDb
from atomate.utils.utils import get_logger
from pymongo import MongoClient, ReturnDocument
from monty.json import jsanitize
from monty.serialization import loadfn
from pymatgen import Structure
from pymatgen import Composition
from pymatgen.transformations.standard_transformations import DiscretizeOccupanciesTransformation

__author__ = 'Matthew Horton'
__credits__ = 'Anubhav Jain, Kiran Mathew, Josua Vieten'
__email__ = 'mkhorton@lbl.gov'
# based on FeffCalcDb

logger = get_logger(__name__)


class McsqsCalcDb(CalcDb):

    def __init__(self, host="localhost", port=27017,
                 database="SQS", collection="tasks",
                 user=None, password=None):
        super(McsqsCalcDb, self).__init__(host, port, database, collection, user, password)

    def mmdb(self, db_file, admin=True):
        """
        Return fixed MMDB 

        Returns:
            MMDb object
        """
        host="localhost"
        port=port=27017
        database="SQS"
        collection="tasks"
        user=None
        password=None

        return host, port, database, collection, user, password


    def build_indexes(self, indexes=None, background=True):
        _indexes = indexes if indexes else ['anonymous_formula']
        for i in _indexes:
            self.collection.create_index(i, background=background)

    def reset(self):
        self.collection.delete_many({})
        self.build_indexes()


    def duplicate_checker(self, disordered_struct, db_file=None, max_denominator=None):
        """
        Check if the disordered structure is already present in the database
        Returns 'True' if a duplicate is found.

        Args:
            disordered_struct: disordered input structure
            db_file: database file
            max_denominator: maximum denominator the occupancies are rounded to

        """

        if db_file:

            # apply the same transformations as in core.py
            if max_denominator:
                trans = DiscretizeOccupanciesTransformation(max_denominator,fix_denominator=True,tol=2)
                disordered_struct = trans.apply_transformation(disordered_struct)
            disordered_struct.scale_lattice(1.0)
            de = disordered_struct

            # look for structure in database
            dbstruct = self.collection.find()

            nextone = 0
            duplicate = None
            while nextone == 0:
                try:
                    record = dbstruct.next()
                    disstruct = record['disordered']
                    struchere = Structure.from_dict(disstruct)
                    if struchere == de: duplicate = struchere
                except:
                    nextone = 1

            if duplicate is None:
                return False

            else:
                return True
   
        else:
            return False


    def inserttask(self, d, update_duplicates=True):
        """
        Insert the task document ot the database collection.

        Args:
            d (dict): task document
            update_duplicates (bool): whether to update the duplicates
        """
        result = self.collection.find_one({"disordered": d["disordered"]}, ["dir_name", "task_id"])
        if result is None or update_duplicates:
            d["last_updated"] = datetime.datetime.today()
            if result is None:
                if ("task_id" not in d) or (not d["task_id"]):
                    d["task_id"] = self.db.counter.find_one_and_update(
                        {"_id": "taskid"}, {"$inc": {"c": 1}},
                        return_document=ReturnDocument.AFTER)["c"]
                logger.info("Inserting {} with taskid = {}".format(d["anonymous_formula"], d["task_id"]))
            elif update_duplicates:
                d["task_id"] = result["task_id"]
                logger.info("Updating {} with taskid = {}".format(d["anonymous_formula"], d["task_id"]))
            d = jsanitize(d, allow_bson=True)
            self.collection.update_one({"dir_name": d["task_id"]},
                                       {"$set": d}, upsert=True)
            return d["task_id"]
        else:
            logger.info("Skipping duplicate {}".format(d["anonymous_formula"]))
            return None

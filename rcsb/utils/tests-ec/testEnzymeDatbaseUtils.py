##
# File:    EnzymeDatabaseUtilsTests.py
# Author:  J. Westbrook
# Date:    3-Feb-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for various utilities for extracting data from Enzyme Database
export data files.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.ec.EnzymeDatabaseUtils import EnzymeDatabaseUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s')
logger = logging.getLogger()


class EnzymeDatabaseUtilsTests(unittest.TestCase):

    def setUp(self):
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), 'rcsb', 'mock-data')
        self.__workPath = os.path.join(HERE, 'test-output')
        #

    def tearDown(self):
        pass

    def testReloadEnzymeDatabase1(self):
        """ Test load from source
        """
        edbu = EnzymeDatabaseUtils(enzymeDirPath=self.__workPath, useCache=False, clearCache=True)
        ecId = '1.2.3.4'
        cl = edbu.getClass(ecId)
        self.assertEqual(cl, 'oxalate oxidase')
        logger.debug("ecId %s class %r" % (ecId, cl))
        linL = edbu.getLineage('1.2.3.4')
        self.assertEqual(len(linL), 4)
        logger.debug("ecId %s lineage (%d) %r" % (ecId, len(linL), linL))

    def testReloadEnzymeDatabase2(self):
        """ Test load from cache
        """
        edbu = EnzymeDatabaseUtils(enzymeDirPath=self.__workPath, useCache=True)
        ecId = '1.2.3.4'
        cl = edbu.getClass(ecId)
        self.assertEqual(cl, 'oxalate oxidase')
        logger.debug("ecId %s class %r" % (ecId, cl))
        linL = edbu.getLineage('1.2.3.4')
        self.assertEqual(len(linL), 4)
        logger.debug("ecId %s lineage (%d) %r" % (ecId, len(linL), linL))


def readEnzymeDatabase():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(EnzymeDatabaseUtilsTests("testReloadEnzymeDatabase1"))
    suiteSelect.addTest(EnzymeDatabaseUtilsTests("testReloadEnzymeDatabase2"))
    return suiteSelect


if __name__ == '__main__':
    if True:
        mySuite = readEnzymeDatabase()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

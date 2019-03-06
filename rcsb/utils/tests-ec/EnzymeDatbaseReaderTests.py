##
# File:    EnzymeDatabaseReaderTests.py
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

from rcsb.utils.ec.EnzymeDatabaseReader import EnzymeDatabaseReader
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s')
logger = logging.getLogger()


class EnzymeDatabaseReaderTests(unittest.TestCase):

    def setUp(self):
        self.__mU = MarshalUtil()
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), 'rcsb', 'mock-data')
        self.__xmlPath = os.path.join(self.__dirPath, 'ec', 'enzyme-data.xml')
        #
        self.__jsonPath = os.path.join(HERE, 'test-output', 'enzyme-lineage-data.json')
        self.__fetchedCopyPath = os.path.join(HERE, 'test-output', 'enzyme-data-copy.xml.gz')
        #

    def tearDown(self):
        pass

    def testFetchDatabase(self):
        edbr = EnzymeDatabaseReader()
        ok = edbr.fetch(self.__fetchedCopyPath)
        self.assertTrue(ok)

    def testReadEnzymeDatabase(self):
        edbr = EnzymeDatabaseReader()
        cD = edbr.read(self.__xmlPath)
        logger.info("Enzyme table data length %d" % len(cD))
        self.assertGreaterEqual(len(cD), 7520)
        ok = self.__mU.doExport(self.__jsonPath, cD, format="json")
        self.assertTrue(ok)


def readEnzymeDatabase():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(EnzymeDatabaseReaderTests("testFetchDatabase"))
    suiteSelect.addTest(EnzymeDatabaseReaderTests("testReadEnzymeDatabase"))

    return suiteSelect


if __name__ == '__main__':
    if True:
        mySuite = readEnzymeDatabase()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

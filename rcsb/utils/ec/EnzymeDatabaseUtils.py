##
# -*- coding: utf-8 -*-
#
# File:    EnzymeDatabaseUtils.py
# Author:  J. Westbrook
# Date:    24-Jan-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Various utilities for extracting data Enzyme database export data files
and returning lineage details.
"""

import copy
import gzip
import logging
import os
import time

from bs4 import BeautifulSoup

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

logger = logging.getLogger(__name__)


class EnzymeDatabaseUtils(object):
    """Various utilities for extracting data Enzyme database export data files
       and returning lineage details.
    """

    def __init__(self, **kwargs):
        #
        urlTarget = kwargs.get('urlTarget', "https://www.enzyme-database.org/downloads/enzyme-data.xml.gz")
        enzymeDirPath = kwargs.get("enzymeDirPath", '.')
        useCache = kwargs.get("useCache", True)
        clearCache = kwargs.get("clearCache", False)
        enzymeDataFileName = kwargs.get('enzymeDataFileName', 'enzyme-data.json')

        self.__debug = False
        #
        self.__mU = MarshalUtil(workPath=enzymeDirPath)
        self.__enzD = self.__reload(urlTarget, enzymeDirPath, enzymeDataFileName=enzymeDataFileName, useCache=useCache, clearCache=clearCache)

    def getClass(self, ecId):
        try:
            return self.__enzD['class'][ecId]
        except Exception:
            pass
        return None

    def getLineage(self, ecId):
        try:
            return self.__enzD['lineage'][ecId]
        except Exception:
            pass
        return None

    def __reload(self, urlTarget, dirPath, enzymeDataFileName, useCache=True, clearCache=False):
        """ Reload input XML database dump file and return data transformed lineage data objects.
'
        Returns:
            dictionary[ec_id] = {'name_list': ... , 'id_list': ... 'depth_list': ... }
        """
        enzD = {}
        #
        fU = FileUtil()
        fn = fU.getFileName(urlTarget)
        xmlFilePath = os.path.join(dirPath, fn)
        enzymeDataPath = os.path.join(dirPath, enzymeDataFileName)
        #
        if clearCache:
            try:
                os.remove(xmlFilePath)
                os.remove(enzymeDataPath)
            except Exception:
                pass
        #
        if useCache and fU.exists(enzymeDataPath):
            enzD = self.__mU.doImport(enzymeDataPath, format='json')
        else:
            if useCache and fU.exists(xmlFilePath):
                logger.debug("Using an existing resource file %s" % xmlFilePath)
                ok = True
            else:
                ok = fU.get(urlTarget, xmlFilePath)
            if ok:
                xrt = self.__parse(xmlFilePath)
                if self.__debug:
                    self.__traverse(xrt, ns="")
                rD = self.__extract(xrt)
                enzD = self.__build(rD)
                ok = self.__mU.doExport(enzymeDataPath, enzD, format='json')
        return enzD

    def __build(self, rD):
        """  Build the list of ancestory classifiers for each leaf classifier.

        Source enzyme database exported data has the following organization:

         "class": [
                    {
                    "id": "1",
                    "class": "1",
                    "subclass": "1",
                    "subsubclass": "0",
                    "heading": "Acting on the CH-OH group of donors",
                    "note": "This subclass contains dehydrogenases that act on primary alcohols,",
                    "last_change": "2006-05-18 10:47:54"
                    }, ...

        "entry": [
                    {
                    "ec_num": "1.1.1.1",
                    "accepted_name": "alcohol dehydrogenase",
                    "reaction": "(1) a primary alcohol + NAD+ = an aldehyde + NADH + H+;;(2) a secondary alcohol + NAD+ = a ketone + NADH + H+",
                    "other_names": "aldehyde reductase; ADH; alcohol dehydrogenase (NAD); aliphatic alcohol dehydrogenase; ethanol dehydrogenase; NAD-dependent alcohol dehydrogenase;",
                    "sys_name": "alcohol:NAD+ oxidoreductase",
                    "comments": "A zinc protein. Acts on primary or secondary alcohols or hemi-acetals with very broad specificity; however the enzyme oxidizes methanol much more poorly than ethanol. ",
                    "links": "BRENDA, EXPASY, GTD, IUBMB, KEGG, PDB, UM-BBD",
                    "class": "1",
                    "subclass": "1",
                    "subsubclass": "1",
                    "serial": "1",
                    "status": null,
                    "diagram": null,
                    "cas_num": "9031-72-5",
                    "glossary": null,
                    "last_change": "2012-02-17 13:19:04",
                    "id": "1"
                    }, ...


        Returns:
            dictionary[ec_id] = {'name_list': ... , 'id_list': ... 'depth_list': ... }

        """
        linD = {}
        cD = {}
        classD = {}
        #
        if 'class' in rD:
            for d in rD['class']:
                clTup = tuple([d[at] for at in ['class', 'subclass', 'subsubclass']])
                val = self.__stripMarkup(d['heading'])
                cD[clTup] = val
                classD['.'.join(clTup)] = val

        for k in sorted(classD):
            logger.debug("%10s %s" % (k, classD[k]))
        #
        if 'entry' in rD:
            for d in rD['entry']:
                if 'serial' not in d:
                    logger.error("Missing serial fieled in %r" % d)
                    continue
                cL = list([d[at] for at in ['class', 'subclass', 'subsubclass', 'serial']])
                #
                tt = cL[:-1]
                sscParent = tuple(tt)
                tt = cL[:-2]
                tt.append("0")
                scParent = tuple(tt)
                tt = cL[:-3]
                tt.extend(["0", "0"])
                cParent = tuple(tt)
                #
                try:
                    idL = ['.'.join(cL[:-3]), '.'.join(cL[:-2]), '.'.join(cL[:-1]), '.'.join(cL)]
                    classVal = cD[cParent] if cParent in cD else ''
                    subClassVal = cD[scParent] if scParent in cD else ''
                    subsubClassVal = cD[sscParent] if sscParent in cD else ''
                    serialValue = d['accepted_name'] if 'accepted_name' in d else ''
                    if serialValue is None or len(serialValue) < 1:
                        serialValue = ''
                        logger.debug("Missing accepted name data for %s" % '.'.join(cL))
                    else:
                        serialValue = self.__stripMarkup(serialValue)
                        classD['.'.join(cL)] = serialValue

                    nmL = [classVal, subClassVal, subsubClassVal, serialValue]
                    depthL = [1, 2, 3, 4]
                    # linD['.'.join(cL)] = {'name_list': nmL, 'id_list': idL, 'depth_list': depthL}
                    linD['.'.join(cL)] = [(depthL[ii], idL[ii], nmL[ii]) for ii in range(4)]
                except Exception as e:
                    logger.error("cL %r cParent %r scParent %r sscParent %r val %r" % (cL, cParent, scParent, sscParent, val))
                    logger.error("d is %r" % (d))
                    logger.error("Failing with %s" % str(e))
        #
        enzD = {"class": classD, 'lineage': linD}
        return enzD

    def __extract(self, xrt):
        """ Extract data from the input document and return a dictionary
            of objects containing rows of dictionaries with attribute naming.

        Args:
            xrt: ElementTree root element

        Returns:
            Extracted data (dict): dictionary organized by category with
                                   XML native data names.
        """
        rD = {}
        for el in xrt.getroot():
            logger.debug("-- Element tag %r name %r" % (el.tag, el.attrib['name']))
            # rD.setdefault(el.tag, []).append(q)
            #
            for ch in el:
                logger.debug("-- --> child element tag %r attrib %r" % (ch.tag, ch.attrib['name']))
                #
                if ch.tag == 'table_data' and ch.attrib['name'] in ['class', 'entry']:
                    dL = self.__getTableData(ch)
                    rD[ch.attrib['name']] = dL

                # for gch in ch:
                #    logger.debug("-- -- --> grand child element tag %r attrib count %r" % (gch.tag, len(gch.attrib)))
                #    for ggch in gch:
                #        logger.debug("-- -- --> grand child element tag %r attrib %r" % (ggch.tag, ggch.attrib['name']))
                    # add parent cardinal attributes
        return rD
        #

    def __getTableData(self, el):
        """ Parse table data sections:

           Example --

             'table_data' {'name': 'cite'}
                'row' {}
                'field' {'name': 'cite_key'}
                'field' {'name': 'ec_num'}
                'field' {'name': 'ref_num'}
                'field' {'name': 'acc_no'}
                'field' {'name': 'last_change'}
                'row' {}
                'field' {'name': 'cite_key'}
                'field' {'name': 'ec_num'}
                'field' {'name': 'ref_num'}
                'field' {'name': 'acc_no'}
                'field' {'name': 'last_change'}
        """
        dL = []
        if el.tag != 'table_data':
            return dL
        for ch in el:
            if ch.tag == 'row':
                d = {}
                for gch in ch:
                    if gch.tag == 'field':
                        nm = gch.attrib['name']
                        val = gch.text
                        d[nm] = val
                        logger.debug(" Field %s val %r" % (nm, val))
                #
                dL.append(copy.copy(d))
        return dL

    def __stripMarkup(self, text):
        try:
            tS = text.replace("&#151;", "-").replace("&amp;#151;", "-").replace(u'â€”', "-")
            self.__bs = BeautifulSoup(tS, features="lxml")
            return self.__bs.get_text()
        except Exception as e:
            logger.exception("Failing for %r with %s" % (text, str(e)))
        return ''

    def __parse(self, filePath):
        """ Parse the input XML data file and return ElementTree root element.
        """
        tree = []
        if filePath[-3:] == '.gz':
            with gzip.open(filePath, mode='rb') as ifh:
                logger.debug('Parsing %s', filePath)
                t = time.time()
                tree = ET.parse(ifh)
                logger.debug('Parsed %s %.2f seconds' % (filePath, time.time() - t))
        else:
            with open(filePath, mode='rb') as ifh:
                logger.debug('Parsing %s', filePath)
                t = time.time()
                tree = ET.parse(ifh)
                logger.debug('Parsed %s in %.2f seconds' % (filePath, time.time() - t))
        return tree

    # -
    def __traverse(self, xrt, ns):
        """ Internal routine to traverse the dom covering/logging all elements and attributes.

        Args:
            xrt (object): ElementTree root element
            ns (str): XML namespace

        """

        for el in xrt.getroot():
            pEl = el.tag.replace(ns, "")
            logger.info("-- %r %r" % (pEl, el.attrib))
            for ch in el:
                chEl = ch.tag.replace(ns, "")
                logger.info("-- -->  %r %r" % (chEl, ch.attrib))
                if ch.text is not None and not len(ch.text):
                    logger.info("-- -->  %s" % ch.text)
                for gch in ch:
                    gchEl = gch.tag.replace(ns, "")
                    logger.info("-- -- -->  %r %r" % (gchEl, gch.attrib))
                    if gch.text is not None and not len(gch.text):
                        logger.info("-- -- -->  %s" % gch.text)
                    for ggch in gch:
                        ggchEl = ggch.tag.replace(ns, "")
                        logger.info("-- -- -- -->  %r %r" % (ggchEl, ggch.attrib))
                        if ggch.text is not None and not len(ggch.text):
                            logger.info("-- -- -- -->  %s" % ggch.text)
                        for gggch in ggch:
                            gggchEl = gggch.tag.replace(ns, "")
                            logger.info("-- -- -- -- -->  %r %r" % (gggchEl, gggch.attrib))
                            if gggch.text is not None and not len(gggch.text):
                                logger.info("-- -- -- -- -->  %s" % gggch.text)

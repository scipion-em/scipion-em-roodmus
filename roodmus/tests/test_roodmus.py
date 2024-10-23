# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     David Herreros (dherreros@cnb.csic.es)
# *
# * National Centre for Biotechnology (CSIC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from pyworkflow.tests import *

from pwem.protocols import ProtImportPdb

from roodmus.protocols import ProtSimulateMicrographs


class TestRoodmusBase(BaseTest):

    def runRoodmus(self, pdbID):
        print("Import atomic model")
        protImportModel = self.newProtocol(ProtImportPdb,
                                           pdbId=pdbID,
                                           objLabel="Refernece model")
        self.launchProtocol(protImportModel)
        self.assertIsNotNone(protImportModel.getStatus(),
                             "There was a problem with the import")

        print("Run Roodmus")
        protSimMic = self.newProtocol(ProtSimulateMicrographs,
                                      topFile=protImportModel.outputPdb,
                                      useGPU=False,
                                      numberOfThreads=4)
        self.launchProtocol(protSimMic)
        self.assertIsNotNone(protSimMic.simMics,
                             "There was a problem with simple initial model protocol")
        self.assertIsNotNone(protSimMic.trueCTFs,
                             "There was a problem with simple initial model protocol")
        self.assertIsNotNone(protSimMic.trueCoords,
                             "There was a problem with simple initial model protocol")


class TestRoodmus(TestRoodmusBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def test_roodmus(self):
        self.runRoodmus("4ake")

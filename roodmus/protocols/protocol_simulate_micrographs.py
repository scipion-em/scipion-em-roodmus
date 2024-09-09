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


import os
from glob import glob

from enum import Enum

from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message
from pyworkflow.object import Set

from pwem.protocols import EMProtocol
from pwem.objects import Micrograph, SetOfMicrographs

from roodmus import Plugin


class outputs(Enum):
    count = SetOfMicrographs


class ProtSimulateMicrographs(EMProtocol):
    """
    Simulation of micrographs with varying conformational variability with Roodmus
    """
    _label = 'simulate micrographs'
    _devStatus = BETA
    _possibleOutputs = outputs

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label=Message.LABEL_INPUT)

        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation.\
                                     Select the one you want to use.")
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")

        form.addParam('topFile', params.PointerParam,
                      pointerClass="AtomStruct",
                      label='Topology file', important=True,
                      help='Atomic model representing the topography (structure) that has been used during the '
                           'molecular dynamics simulation to generate different conformations (no solvent). If no'
                           'trajectory files have been specified, only the conformation represented by this model '
                           'will be used to simulate the micrographs')

        form.addParam('trajFiles', params.PathParam,
                      label='Path to trajectory files', important=True,
                      help='Path to .dcd files obtained from a molecular dynamics simulation. This file is needed '
                           'to simulate particles with varying conformations.')

        form.addParam('numMic', params.IntParam,
                      validators=[params.Positive],
                      default=10,
                      label='Number of micrographs to simulate', important=True)

        form.addParam('numPart', params.IntParam,
                      validators=[params.Positive],
                      default=10,
                      label='Number of particles per micrograph', important=True)

        form.addParam('numConf', params.IntParam,
                      validators=[params.Positive],
                      default=10,
                      condition="trajFiles",
                      label='Number of conformations to sample', important=True)

        form.addParam("pixelSize", params.FloatParam,
                      default=1.0,
                      validators=[params.Positive],
                      label="Micrograph pixel size")

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep(self.sampleConformationsStep)
        self._insertFunctionStep(self.ssimulateMicrographsStep)
        self._insertFunctionStep(self.createOutputStep)

    def sampleConformationsStep(self):
        trajFilesDir = self.trajFiles.get()
        topFile = self.topFile.get().getFileName()
        numConf = self.numConf.get()

        args = (f"--topfile {topFile} --trajfiles_dir {trajFilesDir} --n_conformations {numConf} --tqdm "
                f"--output_dir {self._getExtraPath('simulated_conformations')}")

        program = Plugin.getRoodmusProgram("conformations_sampling")

        self.runJob(program, args)

    def ssimulateMicrographsStep(self):
        numMic = self.numMic.get()
        numPart = self.numPart.get()
        pixelSize = self.pixelSize.get()

        args = (f"--pdb_dir {self._getExtraPath('simulated_conformations')} "
                f"--mrc_dir {self._getExtraPath('simulated_mics')} -n {numMic} -m {numPart} "
                f"--pixel_size {pixelSize} --tqdm "
                f"--nproc {self.numberOfThreads.get()}")

        if self.usesGpu():
            gpuID = [str(elem) for elem in self.getGpuList()][0]
            args += f' --device "gpu" --gpu_id {gpuID}'
        else:
            args += f' --device "cpu"'

        program = Plugin.getRoodmusProgram("run_parakeet")

        self.runJob(program, args)

    def createOutputStep(self):
        pixelSize = self.pixelSize.get()
        outputMics = self._createSetOfMicrographs()

        for micFile in glob(self._getExtraPath(os.path.join('simulated_mics'), "*.mrc")):
            outputMic = Micrograph()
            outputMic.setFileName(micFile)
            outputMic.setSamplingRate(pixelSize)
            outputMics.append(outputMic)

        outputMics.setSamplingRate(pixelSize)

        self._defineOutputs(simMics=outputMics)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        pass

    def _summary(self):
        summary = []

        if self.isFinished():
            numMic = self.simMics.getSize()
            numPart = self.numPart.get()
            pixelSize = self.pixelSize.get()
            numConf = self.numConf.get()
            summary.append(f"A total of {numMic} micrographs have been generated with the following metadata: ")
            summary.append(f"    - Number of particles per micrograph:  {numPart}")
            summary.append(f"    - Number of sampled conformations:  {numConf}")
            summary.append(f"    - Micrograph pixel size: {pixelSize}")
        else:
            summary.append("Simulating micrographs...")

        return summary

    def _methods(self):
        pass

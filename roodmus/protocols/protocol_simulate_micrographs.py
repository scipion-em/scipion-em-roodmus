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
import yaml
from scipy.spatial.transform import Rotation as R

from enum import Enum

from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message, copyFile, getExt, replaceExt
from pyworkflow.object import Set

from pwem.protocols import EMProtocol
from pwem.objects import Micrograph, SetOfMicrographs, CTFModel, Coordinate, Acquisition, SetOfCoordinates

from roodmus import Plugin


class outputs(Enum):
    count = SetOfMicrographs


class ProtSimulateMicrographs(EMProtocol):
    """
    Simulation of micrographs with varying conformational variability with Roodmus
    """
    _label = 'simulate micrographs'
    _devStatus = BETA
    _micModel = ["talos", "krios"]
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
                           'will be used to simulate the micrographs.')

        form.addSection(label="Conformational sampling")

        form.addParam('trajFiles', params.PathParam,
                      label='Path to trajectory files (Optional - see help)', important=False,
                      allowsNull=True,
                      help='Path to .dcd files obtained from a molecular dynamics simulation. This file is needed '
                           'to simulate particles with varying conformations. If not given, the simulated micrographs '
                           'will only include the conformation represented by the topology file.')

        form.addParam('numConf', params.IntParam,
                      validators=[params.Positive],
                      default=10,
                      condition="trajFiles",
                      label='Number of conformations to sample')

        form.addSection(label="Micrograph simulation")

        group = form.addGroup("Micrograph parameters")

        group.addParam('numMic', params.IntParam,
                      validators=[params.Positive],
                      default=10,
                      label='Number of micrographs to simulate', important=True)

        group.addParam('numPart', params.IntParam,
                      validators=[params.Positive],
                      default=10,
                      label='Number of particles per micrograph', important=True)

        group.addParam("pixelSize", params.FloatParam,
                      default=1.0,
                      validators=[params.Positive],
                      label="Micrograph pixel size")

        group.addParam("nX", params.IntParam,
                      default=1000,
                      validators=[params.Positive],
                      label="Micrograph size along X direction")

        group.addParam("nY", params.IntParam,
                      default=1000,
                      validators=[params.Positive],
                      label="Micrograph size along Y direction")

        group.addParam("mag", params.FloatParam,
                       default=50000,
                       experLevel=params.LEVEL_ADVANCED,
                       validators=[params.Positive],
                       label="Magnification rate")

        group.addParam("q0", params.FloatParam,
                       default=0.07,
                       experLevel=params.LEVEL_ADVANCED,
                       validators=[params.Positive],
                       label="Amplitude contrast")

        group = form.addGroup("Micrograph beam")

        group.addParam('dose', params.FloatParam,
                      default=45.0,
                      label='Electron dose (electrons per square angstrom)')

        group = form.addGroup("Ice parameters")

        group.addParam('iceThickness', params.FloatParam,
                      default=500,
                      label='Ice thickness (angstrom)')

        group = form.addGroup("Microscope parameters")

        group.addParam('micModel', params.EnumParam,
                      choices=self._micModel,
                      display=params.EnumParam.DISPLAY_HLIST,
                      default=0,
                      label='Microscope model')

        form.addSection(label="Microscope lens")

        form.addParam('defocusAverage', params.FloatParam,
                      default=-15000,
                      label='Average defocus (angstrom)',
                      help="In CryoEM, this value is negative (underfocus). Positive values (overfocus) are also "
                           "allowed")

        form.addParam('defocusSTD', params.FloatParam,
                      default=5000,
                      label='Defocus standard deviation (angstrom)')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep(self.sampleConformationsStep)
        self._insertFunctionStep(self.simulateMicrographsStep)
        self._insertFunctionStep(self.createOutputStep)

    def sampleConformationsStep(self):
        trajFilesDir = self.trajFiles.get()
        topFile = self.topFile.get().getFileName()
        numConf = self.numConf.get()

        if trajFilesDir:
            args = (f"--topfile {topFile} --trajfiles_dir {trajFilesDir} --n_conformations {numConf} --tqdm "
                    f"--output_dir {self._getExtraPath('simulated_conformations')}")

            program = Plugin.getRoodmusProgram("conformations_sampling")

            self.runJob(program, args)
        else:
            os.mkdir(self._getExtraPath('simulated_conformations'))
            copyFile(topFile, self._getExtraPath(os.path.join('simulated_conformations',
                                                              f"conformation_000000.{getExt(topFile)}")))

    def simulateMicrographsStep(self):
        numMic = self.numMic.get()
        numPart = self.numPart.get()
        pixelSize = self.pixelSize.get()
        iceThickness = self.iceThickness.get()
        nX = self.nX.get()
        nY = self.nY.get()
        centreX = round(0.5 * nX)
        centreY = round(0.5 * nY)
        centreZ = round(0.5 * iceThickness)

        args = (f"--pdb_dir {self._getExtraPath('simulated_conformations')} "
                f"--mrc_dir {self._getExtraPath('simulated_mics')} -n {numMic} -m {numPart} "
                f"--pixel_size {pixelSize} --nx {nX} --ny {nY} --box_x {pixelSize * nX} "
                f"--box_y {pixelSize * nY} --box_z {iceThickness} --centre_x {pixelSize * centreX} "
                f"--centre_y {pixelSize * centreY} --centre_z {centreZ} --cuboid_length_x {pixelSize * nX} "
                f"--cuboid_length_y {pixelSize * nY} --cuboid_length_z {iceThickness} --tqdm "
                f"--nproc {self.numberOfThreads.get()} --electrons_per_angstrom {self.dose.get()} "
                f"--c_10 {self.defocusAverage.get()} --c_10_stddev {self.defocusSTD.get()} ")
                # f"--model {self._micModel[self.micModel.get()]}")  # FIXME: Currently a bug in Roodmus, to be added when fixed

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
        outputCTFs = self._createSetOfCTF()
        outputCoords = self._createSetOfCoordinates(outputMics)
        outputMics.setSamplingRate(pixelSize)

        micId = 1
        for micFile in glob(self._getExtraPath(os.path.join('simulated_mics'), "*.mrc")):
            with open(replaceExt(micFile, "yaml")) as stream:
                yaml_contents = yaml.safe_load(stream)

            # Output 1: Micrographs
            aquisition = Acquisition()
            aquisition.setMagnification(self.mag.get())
            aquisition.setVoltage(yaml_contents["microscope"]["beam"]["energy"])
            aquisition.setDosePerFrame(yaml_contents["microscope"]["beam"]["electrons_per_angstrom"])
            aquisition.setSphericalAberration(yaml_contents["microscope"]["lens"]["c_c"])
            aquisition.setAmplitudeContrast(self.q0.get())
            outputMic = Micrograph()
            outputMic.setFileName(micFile)
            outputMic.setSamplingRate(pixelSize)
            outputMic.setAcquisition(aquisition)
            outputMic.setObjId(micId)
            outputMic.setMicName(f"mic_{micId}")

            # Output 2: CTFs
            ctf = CTFModel()
            ctf.setMicrograph(outputMic)
            ctf.setDefocusU(-yaml_contents["microscope"]["lens"]["c_10"])
            ctf.setDefocusV(-yaml_contents["microscope"]["lens"]["c_10"])
            ctf.setDefocusAngle(yaml_contents["microscope"]["lens"]["phi_12"])
            # outputMic.setCTF(ctf)
            outputCTFs.append(ctf)

            # Output 3: Coordinates
            for pick in yaml_contents["sample"]["molecules"]["local"][0]["instances"]:
                # mat = R.from_euler(angles=pick["orientations"], seq="ZYZ", degrees=False).as_matrix()
                coord = Coordinate()
                coord.setX(int(round(pick["position"][0])))
                coord.setY(int(round(pick["position"][1])))
                coord.setMicrograph(outputMic)
                coord.setMicName(outputMic.getMicName())
                coord.setMicId(outputMic.getObjId())
                outputCoords.append(coord)

            outputMics.append(outputMic)
            outputMics.setAcquisition(aquisition)

            micId += 1

        outputCTFs.setMicrographs(outputMics)
        outputCoords.setMicrographs(outputMics)
        outputCoords.setBoxSize(int(self.nX.get() / 10))

        self._defineOutputs(simMics=outputMics, trueCTFs=outputCTFs, trueCoords=outputCoords)
        self._defineCtfRelation(outputMics, outputCTFs)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        pass

    def _summary(self):
        summary = []

        if self.isFinished():
            numMic = self.simMics.getSize()
            numPart = self.numPart.get() if self.trajFiles.get() else 1
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

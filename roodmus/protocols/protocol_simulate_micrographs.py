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
import numpy as np
import subprocess
from scipy.spatial.transform import Rotation as R

from enum import Enum

from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message, copyFile, getExt, replaceExt
from pyworkflow.object import Set
from pyworkflow.protocol import STEPS_PARALLEL

from pwem.protocols import EMProtocol
from pwem.objects import Micrograph, SetOfMicrographs, CTFModel, Coordinate, Particle, Acquisition, SetOfCoordinates, Transform, String, Boolean
from pwem import ALIGN_PROJ
from pwem.convert import euler_matrix
import pyworkflow.utils as pwutils

from xmipp_metadata.image_handler import ImageHandler

from roodmus import Plugin
from roodmus.utils import normalize_image


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
    stepsExecutionMode = STEPS_PARALLEL

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

        form.addParam('astigmatism', params.NumericListParam,
                      default="10 80",
                      label='The 2-fold astigmatism (angstrom)',
                      help="You can provide here a single value (all micrographs will have the same astigmatism) "
                           "or two values separated by a white spaces (each micrograph will have a random "
                           "astigmatism taken within the range defined by the two numbers provided).")

        form.addSection(label="Output particles")

        form.addParam('boxSize', params.IntParam,
                      default=128,
                      label='Extracted particles box size (px)')

        form.addParam('invertContrast', params.BooleanParam,
                      default=True,
                      label='Invert particle contrast?',
                      help="When set to Yes, particles will be white on a dark background. Otherwise, particles will "
                           "be black in a bright background.")

        form.addParam('doNormalize', params.BooleanParam,
                      default=True,
                      label='Normalize images?',
                      help="Determine whether images are normalized to have zero mean and standard deviation one in the "
                           "background.")

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        deps_preprocessing = self._insertFunctionStep(self.sampleConformationsStep, prerequisites=[])

        deps_simulate = []
        needsGPU = self.usesGpu()
        for idm in range(self.numMic.get()):
            deps_simulate.append(self._insertFunctionStep(self.simulateMicrographsStep,
                                                          idm, prerequisites=deps_preprocessing, needsGPU=needsGPU))

        self._insertFunctionStep(self.createOutputStep, prerequisites=deps_simulate)

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

    def simulateMicrographsStep(self, idm):
        numPart = self.numPart.get()
        pixelSize = self.pixelSize.get()
        iceThickness = self.iceThickness.get()
        nX = self.nX.get()
        nY = self.nY.get()
        centreX = round(0.5 * nX)
        centreY = round(0.5 * nY)
        centreZ = round(0.5 * iceThickness)
        astigmatism = list(map(float, self.astigmatism.get().split(' ')))
        phi_12 = 0.0 if min(astigmatism) == max(astigmatism) == 0.0 else np.random.uniform(0, np.pi)
        astigmatism = np.random.uniform(min(astigmatism), max(astigmatism))

        args = (f"--pdb_dir {self._getExtraPath('simulated_conformations')} "
                f"--mrc_dir {self._getExtraPath(f'simulated_mic_{idm:05}')} -n 1 "
                f"--pixel_size {pixelSize} --nx {nX} --ny {nY} --box_x {pixelSize * nX} "
                f"--box_y {pixelSize * nY} --box_z {iceThickness} --centre_x {pixelSize * centreX} "
                f"--centre_y {pixelSize * centreY} --centre_z {centreZ} --cuboid_length_x {pixelSize * nX} "
                f"--cuboid_length_y {pixelSize * nY} --cuboid_length_z {iceThickness} --tqdm "
                f"--nproc 20 --electrons_per_angstrom {self.dose.get()} "
                f"--c_10 {self.defocusAverage.get()} --c_10_stddev {self.defocusSTD.get()} "
                f"--c_12 {-astigmatism} --phi_12 {phi_12} ")
                # f"--model {self._micModel[self.micModel.get()]}")  # FIXME: Currently a bug in Roodmus, to be added when fixed

        if self.usesGpu():
            args += f' --device gpu --gpu_id  %(GPU)s'
        else:
            args += f' --device cpu'

        program = Plugin.getRoodmusProgram("run_parakeet")
        program_ctf = Plugin.getParakeetProgram("ctf")

        for currNumPart in range(numPart, 1, -10):
            args_with_particles = args +  f' -m {currNumPart}'
            try:
                self.runJob(program, args_with_particles)
                config_file = self._getExtraPath(os.path.join(f'simulated_mic_{idm:05}', f'{0:06}.yaml'))
                ctf_file = self._getExtraPath(os.path.join(f'simulated_mic_{idm:05}', f'{0:06}_ctf.mrc'))
                args_ctf = f'-c {config_file} -o {ctf_file}'
                self.runJob(program_ctf, args_ctf)
                ctf = ImageHandler().read(ctf_file).getData()
                ImageHandler().write(ctf, ctf_file, overwrite=True)
                return
            except subprocess.CalledProcessError as e:
                pwutils.cleanPattern(self._getExtraPath(os.path.join(f'simulated_mic_{idm:05}', "*")))
                if currNumPart - 10 > 0:
                    print(pwutils.yellowStr(f"Could not place the specified number of particles ({numPart}) in "
                                            f"micrograph #{idm}. Retrying Roodmus with {currNumPart - 10} particles"), flush=True)
                else:
                    print(pwutils.redStr(e), flush=True)

    def createOutputStep(self):
        pixelSize = self.pixelSize.get()
        boxSize = self.boxSize.get()
        if boxSize % 2 != 0:
            boxSize += 1
        halfSize = int(round(0.5 * boxSize))
        invertContrast = self.invertContrast.get()
        doNormalize = self.doNormalize.get()
        nX = self.nX.get()
        nY = self.nY.get()
        boxSize = self.boxSize.get()
        outputMics = self._createSetOfMicrographs()
        outputCTFs = self._createSetOfCTF()
        outputCoords = self._createSetOfCoordinates(outputMics)
        outputParticles = self._createSetOfParticles()
        outputMics.setSamplingRate(pixelSize)

        micId = 1
        partId = 1
        particleImgs = []
        stack_file = self._getExtraPath("particle_stack.mrcs")
        for idm in range(self.numMic.get()):
            for micFile in glob(self._getExtraPath(os.path.join(f'simulated_mic_{idm:05}'), "*[!ctf].mrc")):
                with open(replaceExt(micFile, "yaml")) as stream:
                    yaml_contents = yaml.safe_load(stream)

                # Read mic
                micImg = np.squeeze(ImageHandler().read(micFile).getData())

                # Rename file
                if idm > 0:
                    newMicFile = os.path.join(os.path.dirname(micFile), f"{idm:06}.mrc")
                    pwutils.moveFile(micFile, newMicFile)
                else:
                    newMicFile = micFile

                # Output 1: Micrographs
                aquisition = Acquisition()
                aquisition.setMagnification(self.mag.get())
                aquisition.setVoltage(yaml_contents["microscope"]["beam"]["energy"])
                aquisition.setDosePerFrame(yaml_contents["microscope"]["beam"]["electrons_per_angstrom"])
                aquisition.setSphericalAberration(yaml_contents["microscope"]["lens"]["c_c"])
                aquisition.setAmplitudeContrast(self.q0.get())
                outputMic = Micrograph()
                outputMic.setFileName(newMicFile)
                outputMic.setSamplingRate(pixelSize)
                outputMic.setAcquisition(aquisition.clone())
                outputMic.setObjId(micId)
                outputMic.setMicName(f"mic_{micId}")

                # Output 2: CTFs
                ctf = CTFModel()
                ctf.setMicrograph(outputMic)
                ctf.setPsdFile(self._getExtraPath(os.path.join(f'simulated_mic_{idm:05}'), f"{0:06}_ctf.mrc"))
                ctf.setDefocusU(-yaml_contents["microscope"]["lens"]["c_10"] + yaml_contents["microscope"]["lens"]["c_12"])
                ctf.setDefocusV(-yaml_contents["microscope"]["lens"]["c_10"] - yaml_contents["microscope"]["lens"]["c_12"])
                ctf.setDefocusAngle(np.rad2deg(yaml_contents["microscope"]["lens"]["phi_12"]))
                ctf.setPhaseShift(yaml_contents["microscope"]["phase_plate"]["phase_shift"])
                outputMic.setCTF(ctf.clone())
                outputCTFs.append(ctf)

                # Output 3 - 4: Coordinates and Particles
                for pick in yaml_contents["sample"]["molecules"]["local"][0]["instances"]:
                    cx, cy = int(round(pick["position"][0])), int(round(pick["position"][1]))
                    # M = euler_matrix(pick["orientation"][0], pick["orientation"][1], pick["orientation"][2], "szyz")
                    M = np.eye(4)
                    # M[:3, :3] = R.from_euler('ZYZ', pick["orientation"], degrees=False).as_matrix()
                    M[:3, :3] = R.from_rotvec(pick["orientation"]).as_matrix()
                    M[:3, 3] = np.asarray([float(cx) - pick["position"][0], float(cy) - pick["position"][1], 0.0])
                    M = np.linalg.inv(M)
                    tr = Transform()
                    tr.setMatrix(M)
                    coord = Coordinate()
                    coord.setX(cx)
                    coord.setY(cy)
                    coord.setMicrograph(outputMic)
                    coord.setMicName(outputMic.getMicName())
                    coord.setMicId(outputMic.getObjId())
                    outputCoords.append(coord)

                    if ((cx + halfSize < nX) and (cx - halfSize > 0) and
                        (cy + halfSize < nY) and (cy - halfSize > 0)):
                        part = Particle()
                        part.setLocation(partId, stack_file)
                        part.setMicId(outputMic.getObjId())
                        part.setCTF(ctf.clone())
                        part.setTransform(tr.clone())
                        part.setSamplingRate(pixelSize)
                        part.setCoordinate(coord.clone())
                        part.setAcquisition(aquisition.clone())

                        particleImg = micImg[(cy - halfSize):(cy + halfSize), (cx - halfSize):(cx + halfSize)]

                        if invertContrast:
                            particleImg = -1. * particleImg

                        if doNormalize:
                            particleImg = normalize_image(particleImg)

                        particleImgs.append(particleImg)

                        if partId == 1:
                            ImageHandler().write(particleImg[None, :, :], stack_file, overwrite=True)

                        partId += 1

                        outputParticles.append(part)

                outputMics.append(outputMic)
                outputMics.setAcquisition(aquisition.clone())

                micId += 1

        ImageHandler().write(np.stack(particleImgs, axis=0), stack_file, overwrite=True)

        outputCTFs.setMicrographs(outputMics)
        outputCoords.setMicrographs(outputMics)
        outputCoords.setBoxSize(boxSize)
        outputParticles.setSamplingRate(pixelSize)
        outputParticles.setHasCTF(True)
        outputParticles.setAlignmentProj()

        if outputMics.getSize() == 0:
            raise ValueError(pwutils.redStr("No micrographs has been generated by Roodmus. Exiting..."))

        self._defineOutputs(simMics=outputMics, trueCTFs=outputCTFs, trueCoords=outputCoords, trueParticles=outputParticles)
        self._defineCtfRelation(outputMics, outputCTFs)

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

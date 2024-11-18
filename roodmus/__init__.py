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
import pyworkflow.utils as pwutils
import pwem
import subprocess

from roodmus.constants import *

__version__ = "1.0.4"  # plugin version
_logo = "ccpem_logo.png"
_references = ['roodmus2024']

# For the installation
driver_cuda_compatibility = {
    "465": "11.3.0",
    "470": "11.4.0",
    "495": "11.5.0",
    "510": "11.6.0",
    "515": "11.7.0",
    "520": "11.8.0",
    "525": "12.0.0",
    "530": "12.1.0",
    "535": "12.2.0",
    "540": "12.3.0",
    "545": "12.4.0",
    "550": "12.5.0",
    "555": "12.6.0",  # Expected or approximated for newer releases
}


class Plugin(pwem.Plugin):
    _url = "https://github.com/scipion-em/scipion-em-roodmus"
    _supportedVersions = [V1]  # binary version

    @classmethod
    def getEnvActivation(cls):
        return f"conda activate roodmus-{V1}"

    @classmethod
    def getEnviron(cls, gpuID=None):
        """ Setup the environment variables needed to launch my program. """
        environ = pwutils.Environ(os.environ)

        if gpuID is not None:
            environ["CUDA_VISIBLE_DEVICES"] = gpuID

        return environ

    @classmethod
    def getRoodmusProgram(cls, program):
        cmd = '%s %s && roodmus %s' % (cls.getCondaActivationCmd(), cls.getEnvActivation(), program)
        return cmd

    @classmethod
    def getCommand(cls, program, args):
        return cls.getRoodmusProgram(program) + args

    @classmethod
    def defineBinaries(cls, env):

        def getRoodmusInstallationCommands():
            nvidiaNVCC = False
            try:
                nvidiaDriverVer = subprocess.Popen(["nvidia-smi",
                                                    "--query-gpu=driver_version",
                                                    "--format=csv,noheader"],
                                                   env=cls.getEnviron(),
                                                   stdout=subprocess.PIPE
                                                   ).stdout.read().decode('utf-8').split(".")[0]
                nvidiaNVCC = subprocess.run(['nvcc', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            except (ValueError, TypeError, FileNotFoundError):
                print("NVCC not found in your system, installing it in the environment...")

            commands = cls.getCondaActivationCmd() + " "
            if nvidiaNVCC:
                commands += f"conda create -n roodmus-{V1} -c conda-forge fftw python=3.10 -y && "
            else:
                compatible_versions = [cuda for drv, cuda in driver_cuda_compatibility.items() if
                                       drv <= nvidiaDriverVer]
                cudaVersion = max(compatible_versions)
                commands += (
                    f"conda create -n roodmus-{V1} -c conda-forge -c nvidia/label/cuda-{cudaVersion} python=3.10 "
                    f"fftw cuda={cudaVersion} -y && ")
            commands += f"conda activate roodmus-{V1} && "
            commands += "pip install roodmus && pip install openmm && "
            commands += ("git clone https://gitlab.com/ccpem/ccpem-pipeliner.git && "
                         "cd ccpem-pipeliner && git checkout bedbedbe183ad497dbaa82a638f210d316ba9bae && "
                         "pip install -e . && cd .. && ")
            commands += "touch roodmus_installed"
            return commands


        installCmds = [(getRoodmusInstallationCommands(), "roodmus_installed")]
        env.addPackage('roodmus', version=V1,
                       tar='void.tgz',
                       commands=installCmds,
                       default=True)

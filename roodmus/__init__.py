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

from roodmus.constants import *

__version__ = "1.0.2"  # plugin version
_logo = "ccpem_logo.png"
_references = ['roodmus2024']


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
            commands = cls.getCondaActivationCmd() + " "
            commands += f"conda create -n roodmus-{V1} python=3.10 conda-forge::fftw -y && "
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

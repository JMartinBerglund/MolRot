"""
This is the main MolRot module
"""

#!/usr/bin/env python3
from textwrap import dedent
import os
import sys
sys.path.append("/home/martin/Work/GitHub/MolRot/Tomography/")
sys.path.append("/home/martin/Work/GitHub/MolRot/Interferometry/")
sys.path.append("/home/martin/Work/GitHub/MolRot/Utility/")
sys.path.append("/home/martin/Work/GitHub/MolRot/StateCreation/")
sys.path.append("/home/martin/Work/GitHub/MolRot/Vibrations/")
import Interferometry
import Data
import Opt
import Utility
import Tomography
import StateCreation
import Vibrations

print("This is MolRot")
print("Running on python version:")
print(sys.version)
print(sys.version_info)

__all__ = ["Tomography", "Interferometry", "Data", "Opt", "Utility", "StateCreation", "Vibration"]



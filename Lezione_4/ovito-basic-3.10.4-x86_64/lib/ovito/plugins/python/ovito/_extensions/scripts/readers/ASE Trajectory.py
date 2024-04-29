##### ASE trajectory reader #####
#
# [See documentation](manual:file_formats.input.ase_trajectory)

from ovito.data import DataCollection
from ovito.io import FileReaderInterface
from ovito.io.ase import ase_to_ovito
from typing import Callable, Any

# Helper function checking whether the ASE database module is installed in the Python interpreter.
def _load_ase_module():
    try:
        import ase.io
    except ImportError as exc:
        raise ImportError('This file reader requires the ASE Python module. '
            'Please install it first by running "pip install ase" if you are using an external Python interpreter, '
            'or "conda install ase" if you are using an Anaconda environment, '
            'or "ovitos -m pip install --user ase" if you are using the embedded interpreter of OVITO Pro.') from exc
    return ase.io

class ASETrajectoryReader(FileReaderInterface):

    @staticmethod
    def detect(filename: str):
        try:
            aseio = _load_ase_module()
            # This will raise an exception if the file is not a valid ASE trajectory file.
            aseio.read(filename, format="traj")
            return True
        except:
            pass
        try:
            aseio = _load_ase_module()
            # This will raise an exception if the file is not a valid ASE trajectory file.
            aseio.read(filename, format="bundletrajectory")
            return True
        except:
            pass
        return False

    def scan(self, filename: str, register_frame: Callable[..., None]):
        try:
            aseio = _load_ase_module()
            trj = aseio.read(filename, index=':')
            for i in range(len(trj)):
                register_frame()
        except aseio.formats.UnknownFileTypeError as ex:
            raise RuntimeError("Unknown ASE file type") from ex

    def parse(self, data: DataCollection, filename: str, frame_index: int, **kwargs: Any):
        try:
            aseio = _load_ase_module()
            ase_atoms = aseio.read(filename, index=frame_index)
            # Convert the ASE Atoms object to the canonical OVITO representation.
            ase_to_ovito(ase_atoms, data)
        except aseio.formats.UnknownFileTypeError as ex:
            raise RuntimeError("Unknown ASE file type") from ex

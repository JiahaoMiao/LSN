##### ASE database reader #####
#
# [See documentation](manual:file_formats.input.ase_database)

from ovito.data import DataCollection
from ovito.io import FileReaderInterface
from ovito.io.ase import ase_to_ovito
from typing import Callable, Any
import traits.api

# Helper function checking whether the ASE database module is installed in the Python interpreter.
def _load_ase_module():
    try:
        import ase.db
    except ImportError as exc:
        raise ImportError('This file reader requires the ASE Python module. '
            'Please install it first by running "pip install ase" if you are using an external Python interpreter, '
            'or "conda install ase" if you are using an Anaconda environment, '
            'or "ovitos -m pip install --user ase" if you are using the embedded interpreter of OVITO Pro.') from exc
    return ase.db

class ASEDatabaseReader(FileReaderInterface):

    # Optional filter expression selecting only certain structures from the ASE database.
    query_string = traits.api.Str(label='ASE query string')

    @staticmethod
    def detect(filename: str):
        try:
            asedb = _load_ase_module()
            # This will raise an exception if the file is not a valid ASE database file.
            asedb.connect(filename).count()
        except:
            return False
        return True

    def scan(self, filename: str, register_frame: Callable[..., None]):
        # Open the ASE database file selected by the user and
        # step through the rows of the database treating each as a separate trajectory frame.
        asedb = _load_ase_module()
        db = asedb.connect(filename)
        for row in db.select(self.query_string, verbosity=1):
            register_frame(frame_info=row.id, label=row.formula)

    def parse(self, data: DataCollection, filename: str, frame_info: int, **kwargs: Any):
        asedb = _load_ase_module()
        db = asedb.connect(filename)
        # Load the db row corresponding to the current frame number.
        # Its db ID was stored in the 'parser_data' field by the scan() method.
        ase_atoms = db.get_atoms(frame_info, add_additional_information=True)
        # Convert the ASE Atoms object to the canonical OVITO representation.
        ase_to_ovito(ase_atoms, data)

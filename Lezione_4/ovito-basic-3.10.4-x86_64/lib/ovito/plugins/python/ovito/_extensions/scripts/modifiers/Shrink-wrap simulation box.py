##### Shrink-wrap simulation box #####
#
# This modifier resets the simulation cell to the tightly fit bounding box enclosing all particle coordinates.
# [See documentation...](manual:modifiers.shrink_wrap_box).

from ovito.data import DataCollection
import numpy

def modify(frame: int, data: DataCollection):

    # There's nothing we can do if there are no input particles.
    if not data.particles or data.particles.count == 0:
        return

    # Compute min/max range of particle coordinates.
    coords_min = numpy.amin(data.particles.positions, axis=0)
    coords_max = numpy.amax(data.particles.positions, axis=0)

    # Build the new 3x4 cell matrix:
    #   (x_max-x_min  0            0            x_min)
    #   (0            y_max-y_min  0            y_min)
    #   (0            0            z_max-z_min  z_min)
    matrix = numpy.empty((3,4))
    matrix[:,:3] = numpy.diag(coords_max - coords_min)
    matrix[:, 3] = coords_min

    # Assign the cell matrix - or create whole new SimulationCell object in
    # the DataCollection if there isn't one already.
    data.create_cell(matrix, (False, False, False))
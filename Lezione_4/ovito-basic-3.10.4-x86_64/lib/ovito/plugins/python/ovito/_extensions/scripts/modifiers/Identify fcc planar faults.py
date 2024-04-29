##### Identify fcc planar faults #####
#
# This modifier function attributes hcp-like atoms in fcc crystals to different planar defect types such as stacking faults and twin boundaries.
# [See documentation...](manual:modifiers.identify_fcc_planar_faults).

from ovito.data import DataCollection, ElementType, DataTable, PTMNeighborFinder
from ovito.modifiers import PolyhedralTemplateMatchingModifier
from ovito.pipeline import ModifierInterface
from ovito.traits import Color
import ovito
import numpy as np
from enum import IntEnum

class IdentifyFCCPlanarFaultsModifier(ModifierInterface):
    """
    Base: :py:class:`ovito.pipeline.ModifierInterface`

    This :ref:`Python-based modifier <writing_custom_modifiers>` identifies different planar defect types in face-centered cubic (fcc) crystals
    such as stacking faults and coherent twin boundaries. All of these planar defects have in common that they are
    made of hcp-like atoms arranged in {111} planes of the fcc crystal.
    See also the corresponding :ref:`user manual page <manual:modifiers.identify_fcc_planar_faults>` for more information.

    The identification algorithm relies on intermediate results provided by the :py:class:`PolyhedralTemplateMatchingModifier`, which
    must be inserted into the data pipeline first to identify all hcp-like defect atoms in the fcc crystal:

    .. literalinclude:: ../example_snippets/identify_fcc_planar_faults_modifier.py
       :lines: 6-8

    The modifier subsequently analyzes the neighborhood of each hcp-like atom to identify which kind
    of planar fault it is part of. Each atom in the input system is attributed to one of following classes:

      * 0 = Non-hcp structure (e.g. perfect fcc)
      * 1 = Indeterminate hcp-like (e.g. isolated hcp-like atoms, not a planar defect)
      * 2 = Intrinsic stacking fault (two hcp-like layers)
      * 3 = Coherent twin boundary (one hcp-like layer)
      * 4 = Multilayer stacking fault (three or more hcp-like layers)

    The modifier writes the corresponding numeric values to the ``Planar Fault Type`` output particle property.
    Furthermore, the modifier produces a :py:class:`~ovito.data.DataTable` which lists the total atom count and
    estimated total area for each planar defect type. This information can be accessed as follows:

    .. literalinclude:: ../example_snippets/identify_fcc_planar_faults_modifier.py
       :lines: 14-19

    .. versionadded:: 3.8.0
    """

    class Types(IntEnum):
        "The possible classifications for hcp-like defect atoms."
        NONHCP = 0
        "Non-hcp structure"
        OTHER = 1
        "Indeterminate hcp-like"
        ISF = 2
        "Intrinsic stacking fault"
        TWIN = 3
        "Coherent twin boundary"
        MULTI = 4
        "Multilayer stacking fault"

    # Color value associated with non-hcp atoms. This parameter is not really used, because the modifier
    # does not overwrite the existing colors of non-hcp atoms.
    color_nonhcp = Color(label='Non-hcp structure (0)', default=(1.0,1.0,1.0))

    # Define modifier parameters controlling the colors of different hcp-like atoms:
    color_other = Color(label='Indeterminate hcp-like (1)', default=(1.0,0.4,0.4))
    """
    Color value to be assigned to all hcp-like atoms that are not part of an identified planar fault (:py:attr:`Types.OTHER`).

    :Default: ``(1.0, 0.4, 0.4)``
    """

    color_isf = Color(label='Intrinsic stacking fault (2)', default=(0.9,0.7,0.2))
    """
    Color value to be assigned to all hcp-like atoms that are part of an intrinsic stacking fault (:py:attr:`Types.ISF`).

    :Default: ``(0.9, 0.7, 0.2)``
    """

    color_twin = Color(label='Coherent twin boundary (3)', default=(0.4,0.4,1.0))
    """
    Color value to be assigned to all hcp-like atoms that are part of a coherent twin boundary (:py:attr:`Types.TWIN`).

    :Default: ``(0.4, 0.4, 1.0)``
    """

    color_multi = Color(label='Multilayer stacking fault (4)', default=(0.0,1.0,0.9))
    """
    Color value to be assigned to all hcp-like atoms that are part of a multi-layer stacking fault (:py:attr:`Types.MULTI`).

    :Default: ``(0.0, 1.0, 0.9)``
    """

    def modify(self, data: DataCollection, **kwargs):

        if data.particles is None:
            raise RuntimeError('This modifier operates on particles but there are no particles in the current dataset.')

        # First, make sure that a PTM modifier is present in the upstream pipeline and has already identified FCC and HCP atoms.
        # If this is not the case, ask user to insert the PTM modifier into the pipeline before this Python modifier.
        if not 'Structure Type' in data.particles or not 'Correspondences' in data.particles:
            raise RuntimeError('This modifier relies on information produced by the PTM algorithm. '
                'Please insert the Polyhedral Template Matching modifier into the pipeline before this modifier and let it '
                'identify FCC and HCP atoms at least. Make sure the output of lattice orientations and interatomic distances '
                'is enabled in the PTM modifier.')
        structure_types = data.particles.structure_types[...]

        # Zero-based indices of all HCP atoms.
        hcp_indices = np.flatnonzero(structure_types == PolyhedralTemplateMatchingModifier.Type.HCP)

        # Total number of HCP atoms.
        num_hcp = len(hcp_indices)

        # Initialize the PTM neighbor finder facility.
        yield 'Building hcp atom neighbor lists'
        neigh_finder = PTMNeighborFinder(data)

        # Precompute the neighbor lists of all HCP atoms.
        if True:
            # Use vectorized PTMNeighborFinder.find_all() method, which returns an (N,12) array.
            # It contains the list of 12 neighbors of each HCP atom (in standard PTM order).
            hcp_neighbors = neigh_finder.find_all(PolyhedralTemplateMatchingModifier.Type.HCP)
            assert hcp_neighbors.shape == (num_hcp, 12)

            # For all neighbors of HCP atoms that are also HCP, encode the neighbors as a non-negative indices
            # into the hcp_indices array, which contains all HCP atoms.
            neighbors_which_are_hcp = (structure_types[hcp_neighbors] == PolyhedralTemplateMatchingModifier.Type.HCP)
            hcp_neighbors[neighbors_which_are_hcp] = np.searchsorted(hcp_indices, hcp_neighbors[neighbors_which_are_hcp])

            # For all neighbors that are not HCP, encode the indices into the global particles list (i.e., not into the hcp_indices array).
            # To differentiate them from the HCP neighbors above, use negative values.
            neighbors_which_are_not_hcp = np.logical_not(neighbors_which_are_hcp)
            hcp_neighbors[neighbors_which_are_not_hcp] = -hcp_neighbors[neighbors_which_are_not_hcp] - 1
        else:
            # Old implementation of the algorithm (no longer used - just for reference), which was based on a
            # for-loop (slow!) calling the PTMNeighborFinder.find() method for each HCP atom.
            hcp_neighbors = np.empty((num_hcp, 12), dtype=np.int64)
            for i, aindex in enumerate(hcp_indices):
                # Show progress of the algorithm in the OVITO Pro status bar.
                yield i / num_hcp
                # Visit each of the 12 neighbors of the current HCP atom in fixed order.
                for j, neigh in enumerate(neigh_finder.find(aindex)):
                    # Is this neighbor an HCP atom too?
                    if structure_types[neigh.index] == PolyhedralTemplateMatchingModifier.Type.HCP:
                        # If the neighbor of an HCP atom is also an HCP atom, the neighbor atom's index will be
                        # encoded as a non-negative index into the hcp_indices array.
                        hcp_neighbors[i,j] = np.searchsorted(hcp_indices, neigh.index)
                    else:
                        # If the neighbor is not an HCP atom, the atom's index into the global particles list
                        # will be encoded as a negative value.
                        hcp_neighbors[i,j] = -neigh.index - 1 # Encode the non-HCP neighbor's global index as a negative value

        if num_hcp > 0:
            # Initialize the array 'layer_dir' (12 integers) indicating for each of the 12 neighbors of any HCP atom,
            # whether the neighbor is below, within, or above the HCP basal plane (values -1,0,+1).
            # This reflects the hard-coded order in which PTMNeighborFinder visits the 12 neighbors of any HCP atom.
            idealvecs = [neigh.idealvec for neigh in neigh_finder.find(hcp_indices[0])]
            layer_dir = np.array([0 if (v[2] == 0.0) else (1 if (v[2]) > 0.0 else -1) for v in idealvecs])

            # The 6 indices of neighbors within and off the basal plane, as separate arrays.
            basal_neighbors = np.flatnonzero(layer_dir == 0)
            outofplane_neighbors = np.flatnonzero(layer_dir != 0)

        # Allocate output per-particle array, which will store the identified type of each HCP atom.
        fault_types = np.zeros_like(structure_types)

        # ElementType objects for the possible hcp atom types, which will be associated with the "Planar Fault Type" output particle property.
        nonhcp_type = ElementType(id=IdentifyFCCPlanarFaultsModifier.Types.NONHCP, name='Non-hcp',                   color=self.color_nonhcp)
        other_type  = ElementType(id=IdentifyFCCPlanarFaultsModifier.Types.OTHER,  name='Indeterminate hcp-like',    color=self.color_other)
        isf_type    = ElementType(id=IdentifyFCCPlanarFaultsModifier.Types.ISF,    name='Intrinsic stacking fault',  color=self.color_isf)
        twin_type   = ElementType(id=IdentifyFCCPlanarFaultsModifier.Types.TWIN,   name='Coherent twin boundary',    color=self.color_twin)
        multi_type  = ElementType(id=IdentifyFCCPlanarFaultsModifier.Types.MULTI,  name='Multilayer stacking fault', color=self.color_multi)

        # Helper function which determines for two neighboring HCP atoms 'a' and 'b' whether they are properly stacked.
        # Properly stacked means their basal planes are parallel. They are not properly stacked if there is a third atom,
        # which is a common neighbor of 'a' and 'b' and which is in the basal planes of both 'a' and 'b'.
        def are_stacked(a: int, b: int):
            # Iterate over the basal plane neighbors of both 'a' and 'b'.
            for i in basal_neighbors:
                for j in basal_neighbors:
                    if hcp_neighbors[a,i] == hcp_neighbors[b,j]:
                        # The two HCP atoms 'a' and 'b' share a common neighbor, which is both basal planes.
                        return False
            return True

        # First pass over each HCP atom.
        yield 'Identifying planar faults (first pass)'
        for i, aindex in enumerate(hcp_indices):
            # Display progress of algorithm in the OVITO status bar.
            yield i / num_hcp

            n_basal = 0 # Number of HCP neighbors in the same basal plane
            n_positive = 0 # Number of HCP neighbors on one side of the basal plane
            n_negative = 0 # Number of HCP neighbors on opposite side
            n_fcc_positive = 0 # Number of FCC neighbors on one side of the basal plane
            n_fcc_negative = 0 # Number of FCC neighbors on opposite side of the basal plane

            # Visit each of the 12 neighbors of the HCP atom in order.
            for j in range(12):

                # Is the current neighbor an HCP atom as well?
                if hcp_neighbors[i,j] >= 0:

                    # Classify the neighbor as being either in the basal plane or out-of-plane.
                    if layer_dir[j] == 0: n_basal += 1
                    elif are_stacked(i, hcp_neighbors[i,j]):
                        if layer_dir[j] == 1: n_positive += 1
                        else: n_negative += 1

                elif layer_dir[j] != 0:
                    # Check if the HCP atom has FCC neighbors on either side of the basal plane.
                    neighbor_type = structure_types[-hcp_neighbors[i,j] - 1]
                    if neighbor_type == PolyhedralTemplateMatchingModifier.Type.FCC:
                        if layer_dir[j] > 0: n_fcc_positive += 1
                        else: n_fcc_negative += 1

            # Initial classification of HCP atom:

            # Is it an intrinsic stacking fault (two parallel HCP layers, atom has at least one out-of-plane neighbor in one direction)?
            if (n_positive != 0 and n_negative == 0) or (n_positive == 0 and n_negative != 0):
                fault_types[aindex] = isf_type.id
            # Is it a coherent twin boundary (single HCP layer, atom has no out-of-plane HCP neighbors but at least one in-plane neighbor)?
            elif n_basal >= 1 and n_positive == 0 and n_negative == 0 and n_fcc_positive != 0 and n_fcc_negative != 0:
                fault_types[aindex] = twin_type.id
            # Is it a multi-layered stacking fault (three or more HCP layers, atom has out-of-plane HCP neighbors on both sides)?
            elif n_positive != 0 and n_negative != 0:
                fault_types[aindex] = multi_type.id
            # Otherwise, it must be an isolated HCP atom (undetermined planar fault type)
            else:
                fault_types[aindex] = other_type.id

        # Second pass over all HCP atoms.
        yield 'Identifying planar faults (second pass)'
        for i, aindex in enumerate(hcp_indices):
            # Display progress of algorithm in the OVITO status bar.
            yield i / num_hcp

            # Reconsider the classification of coherent twin and undetermined HCP atom.
            if fault_types[aindex] == twin_type.id or fault_types[aindex] == other_type.id:
                n_isf_neighbors = 0 # Counts the number of neighbors in the basal plane which are intrinsic stacking fault (ISF) atoms
                n_twin_neighbors = 0 # Counts the number of neighbors in the basal plane which are coherent twin boundary (CTB) atoms
                # Visit the 6 neighbors in the basal plane.
                for j in basal_neighbors:
                    neighbor_index = hcp_neighbors[i,j]
                    # Is the current neighbor an HCP atom?
                    if neighbor_index >= 0:
                        # Check the planar fault type of the neighbor atom.
                        neighbor_index = hcp_indices[neighbor_index]
                        if fault_types[neighbor_index] == isf_type.id:
                            n_isf_neighbors += 1
                        elif fault_types[neighbor_index] == twin_type.id:
                            n_twin_neighbors += 1
                # If the TB atom is surrounded by ISF atoms only, turn it into an ISF atom as well.
                if n_isf_neighbors != 0 and n_twin_neighbors == 0:
                    fault_types[aindex] = isf_type.id
                # If the Other atom is surrounded by TB atoms only, turn it into a TB atom as well.
                elif n_isf_neighbors == 0 and n_twin_neighbors != 0:
                    fault_types[aindex] = twin_type.id

            elif fault_types[aindex] == multi_type.id:
                # Visit the 6 out-of-plane neighbors.
                for j in outofplane_neighbors:
                    neighbor_index = hcp_neighbors[i,j]
                    # Is the current neighbor an HCP atom of type ISF?
                    if neighbor_index >= 0 and fault_types[hcp_indices[neighbor_index]] == isf_type.id:
                        # Turn the neighbor into a multi-layered fault atom.
                        fault_types[hcp_indices[neighbor_index]] = multi_type.id

        # Create the output particle property 'Planar Fault Type'
        property = data.particles_.create_property('Planar Fault Type', data=fault_types)
        property.types.append(nonhcp_type)
        property.types.append(other_type)
        property.types.append(isf_type)
        property.types.append(twin_type)
        property.types.append(multi_type)

        # Give HCP atoms new coloring based on the identified planar fault types.
        data.particles_['Color_'][fault_types == other_type.id] = other_type.color
        data.particles_['Color_'][fault_types == isf_type.id]   = isf_type.color
        data.particles_['Color_'][fault_types == twin_type.id]  = twin_type.color
        data.particles_['Color_'][fault_types == multi_type.id] = multi_type.color

        # Output identification statistics as a bar chart.
        table = DataTable(plot_mode=DataTable.PlotMode.BarChart, identifier='planar_faults', title='Planar fault atoms', count=5)
        table.x = table.create_property('Planar Fault Type', dtype=int, components=1)
        table.y = table.create_property('Atom Count', dtype=np.int64, components=1)
        fractions = table.create_property('Atom Fraction', dtype=float, components=1)
        for t in (other_type, isf_type, twin_type, multi_type):
            table.x[t.id] = t.id
            table.x.types.append(t)
            table.y[t.id] = np.count_nonzero(fault_types == t.id)
            fractions[t.id] = table.y[t.id] / max(num_hcp, 1)
        data.objects.append(table)

        # Compute estimated areas of the ISF and TB planar fault types.
        if 'Interatomic Distance' in data.particles:
            # Nearest neighbor distances of HCP atoms:
            interatomic_distance = data.particles['Interatomic Distance']
            areas = table.create_property('Estimated Area', dtype=float, components=1)
            # Area projected on {111} plane of one hcp atom is sqrt(3) / 2 * a^2.
            areas[isf_type.id] = np.sqrt(3) / 4 * np.sum(interatomic_distance[fault_types == isf_type.id]**2)
            areas[twin_type.id] = np.sqrt(3) / 2 * np.sum(interatomic_distance[fault_types == twin_type.id]**2)

        if ovito.gui_mode:
            print(f"Analyzed {num_hcp} hcp atoms.")

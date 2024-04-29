import importlib
import ovito.modifiers

# Inject Python script extensions into the ovito.modifiers package.
# Note: Using importlib.import_module() here to load them, because human-readable Python filenames contain whitespace.

ovito.modifiers.IdentifyFCCPlanarFaultsModifier = importlib.import_module(".Identify fcc planar faults", __name__).IdentifyFCCPlanarFaultsModifier
ovito.modifiers.__all__ += ['IdentifyFCCPlanarFaultsModifier']

ovito.modifiers.RenderLAMMPSRegionsModifier = importlib.import_module(".Render LAMMPS regions", __name__).RenderLAMMPSRegionsModifier
ovito.modifiers.__all__ += ['RenderLAMMPSRegionsModifier']

# The following Python modifiers are still implemented as standalone modify() functions, not as modern ModifierInterface classes.
# We import the modify() function into the ovito.modifiers module to make it discoverable by the entry points mechanism (which may be confused by the whitespace in the original module name).
ovito.modifiers.CalculateLocalEntropyFunction = importlib.import_module(".Calculate local entropy", __name__).modify
ovito.modifiers.CalculateLocalEntropyFunction.__name__ = 'CalculateLocalEntropyFunction' # Override function name to replace "modify" with the public name. This is needed for correct Python codegen.
ovito.modifiers.__all__ += ['CalculateLocalEntropyFunction']

ovito.modifiers.ShrinkWrapSimulationBoxFunction = importlib.import_module(".Shrink-wrap simulation box", __name__).modify
ovito.modifiers.ShrinkWrapSimulationBoxFunction.__name__ = 'ShrinkWrapSimulationBoxFunction' # Override function name to replace "modify" with the public name. This is needed for correct Python codegen.
ovito.modifiers.__all__ += ['ShrinkWrapSimulationBoxFunction']

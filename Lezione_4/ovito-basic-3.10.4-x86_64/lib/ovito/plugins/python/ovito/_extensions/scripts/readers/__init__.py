import importlib

# Note: Using importlib.import_module() to import modules, because human-readable Python filenames contain whitespace.
ASEDatabaseReader = importlib.import_module(".ASE Database", __name__).ASEDatabaseReader
ASETrajectoryReader = importlib.import_module(".ASE Trajectory", __name__).ASETrajectoryReader

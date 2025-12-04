from pylib.pysim import Pysim


def reproduce():
    """
    used to reproduce simulation results.
    runs simulation from a directory that already contains:
    - config.json
    - config.json["bead_type_filenames"]
    """

    sim = Pysim.from_directory(".")
    sim.run("reproduce-simulation")


if __name__ == "__main__":
    reproduce()

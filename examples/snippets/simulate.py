from pylib import analysis
from pylib.pysim import Pysim
from pylib.utils import cd

sim = Pysim.from_directory("1024/iteration0")
sim.config["load_configuration"] = False
sim.config["grid_size"] = 500

root = "testing-grid500"
sim.set_root(root, mkdir=True)
sim.setup_needed = True
sim.run("production_out")

with cd(root):
    analysis.main()

from scribe import analysis
from scribe.scribe_sim import ScribeSim
from scribe.utils import cd

sim = ScribeSim.from_directory("1024/iteration0")
sim.config["load_configuration"] = False
sim.config["grid_size"] = 500

root = "testing-grid500"
sim.set_root(root, mkdir=True)
sim.setup_needed = True
sim.run("production_out")

with cd(root):
    analysis.main()

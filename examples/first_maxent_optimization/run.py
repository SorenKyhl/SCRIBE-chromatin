from pylib import default, hic, pipeline

"""
executes simplest maxent optimization using pipeline class and library defaults.
"""

config = default.config
params = default.params
contactmap = hic.load_hic(config["nbeads"])
pipe = pipeline.Pipeline(
    name="my_first_maxent_optimization", gthic=contactmap, config=config, params=params
)
pipe.fit()

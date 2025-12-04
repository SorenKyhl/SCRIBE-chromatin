from pylib.maxent import Maxent

me = Maxent.load_state("backup.pickle")
me.set_root("newroot")
me.fit()

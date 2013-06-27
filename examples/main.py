from kt_simul.core import simul_spindle as sim
from kt_simul.io import SimuIO
from kt_simul.draw import Drawer

meta = sim.Metaphase(verbose = True)
meta.simul()

d = Drawer(meta)
d.show_all(fname = "trajs.pdf")
d.show_one(fname = "one.png")

io = SimuIO(meta)
io.save("results.xml", "data.npy")

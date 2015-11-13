from ..mecabio import Model
from ..mecabio import viscous
from ..mecabio import dampedspring
from ..mecabio import linear_fv


class SpindleModel(Model):

    def __init__(self, spindle):
        self.params = spindle.params
        Model.__init__(self, spindle, dt=1) #, self.params['dt'])
        self.prng = spindle.prng
        self.spindle = spindle
        self.duration = self.params['span']
        self.dt = self.params['dt']
        self.num_steps = int(self.duration / self.dt)
        self.time_invariantA()
        self.anaphase = False

    def time_invariantA(self):

        viscous(self, self.spindle.spbL, self.params['mus'])
        viscous(self, self.spindle.spbR, self.params['mus'])

        # Here friction with nucleoplasm is shared btw all objects

        for ch in self.spindle.chromosomes:

            viscous(self, ch.cen_A, self.params['muco'])
            viscous(self, ch.cen_B, self.params['muco'])

            for ps in ch.cen_A.plugsites:
                viscous(self, ps, self.params['muco'])
            for ps in ch.cen_B.plugsites:
                viscous(self, ps, self.params['muco'])

        self.A0mat = self.Amat.copy()

    def update_AB(self):

        self.Amat[:] = self.A0mat
        self.Bvect *= 0

        mz = self.spindle.links[(self.spindle.spbL.idx, self.spindle.spbR.idx)]
        linear_fv(self, mz, self.params['Fmz'], self.params['Vmz'], gamma=1)

        for ch in self.spindle.chromosomes:

            chromatin = self.spindle.links[(ch.cen_A.idx, ch.cen_B.idx)]
            dampedspring(self, chromatin,
                         self.params['muc'],
                         self.params['kappa_c'],
                         self.params['d0'])

            for ps in ch.cen_A.plugsites:
                self.plugsite_forces(ps)
            for ps in ch.cen_B.plugsites:
                self.plugsite_forces(ps)

    def plugsite_forces(self, ps):

        kt = self.spindle.links[(ps.centromere.idx, ps.idx)]
        dampedspring(self, kt,
                     self.params['muk'],
                     self.params['kappa_k'], 0.01)

        if ps.plug_state == 0:
            return

        if ps.plug_state == -1:
            ktMT = self.spindle.links[(self.spindle.spbL.idx, ps.idx)]
        elif ps.plug_state == 1:
            ktMT = self.spindle.links[(self.spindle.spbR.idx, ps.idx)]

        linear_fv(self, ktMT, 1 * ps.calc_ldep(), 1, -1)

    def one_step(self, step):

        if not self.anaphase:
            for plugsite in self.spindle.all_plugsites():
                plugsite.plug_unplug()

        self.update_AB()
        self.solve()
        self.spindle.register_history(step)

        if step == self.num_steps - 1:
            self.spindle.end_history()

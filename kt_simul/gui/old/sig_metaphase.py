from PySide import QtGui, QtCore

from kt_simul.core.simul_spindle import Metaphase


class SigMetaphase(QtGui.QWidget, Metaphase):

    """
    The aim of this hybrid is to retrieve signals from the
    simulation while it's running.
    overrides _one_step method of the Metaphase class
    """

    def __init__(self, paramtree, measuretree,
                 initial_plug=None, parent=None):

        Metaphase.__init__(self, paramtree, measuretree,
                           initial_plug=initial_plug)
        QtGui.QWidget.__init__(self, parent)

        self.date = 0

    def _one_step(self):

        if not self.KD.anaphase:
            self.emit(QtCore.SIGNAL('inMetaphase'))

        self.KD.one_step(self.date)

        nb_mero = self._mero_checkpoint()
        if nb_mero > 0:
            self.emit(QtCore.SIGNAL('meroCheckPoint'), nb_mero)

        self.emit(QtCore.SIGNAL('plugCheckPoint'), self._plug_checkpoint())
        self.date += 1
        self.emit(QtCore.SIGNAL('stepDone'), self.date)
        self.emit(QtCore.SIGNAL('stepDone_nop'))

    def sig_simul(self):

        self.simul()
        self.emit(QtCore.SIGNAL('simulDone(bool)'), self.report)
        self.simulDone = True

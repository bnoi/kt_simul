# -*- coding: utf-8 -*-

from PySide import QtGui

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

class MyMplCanvas(FigureCanvas):
    """
    Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.).
    """

    def __init__(self, span=800., parent=None,
                    width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure(span)

        #
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self, span):
        self.axes.axis([0, span , -10, 10 ])

        self.axes.set_xlabel('Time (seconds)', fontsize = 12)
        self.axes.set_ylabel(u'Distance from center (µm)', fontsize = 12)

    def update_figure(self, mt, n = None):

        '''
        Plot the different trajectories
        '''

        self.axes.clear()
        self.axes.hold(True)
        if n == None:
            mt.show_trajs(self.axes)
        else:
            mt.show_one(n = n, fig = self.fig)

        self.draw()

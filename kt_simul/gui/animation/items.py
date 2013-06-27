#-*- coding: utf-8 -*-

from PySide import QtCore, QtGui

SPB_WIDTH = 0.5
SPB_HEIGHT = 0.7

PLUGSITE_OFFSET = 0.05  # Vertical distance between attachment sites
CH_OFFSET = 0.4  # Vertical distance between chromosomes

PLUGSITE_WIDTH = 0.1
PLUGSITE_HEIGHT = 0.1

SPB_COLOR = QtGui.QColor(255, 0, 0)  # Red

CH_COLOR = QtGui.QColor(0, 100, 100, alpha=200)  # Green
ACTIVE_SAC_COLOR = QtGui.QColor(100, 10, 100, alpha=200)  # Purple

GOOD_PLUGSITE_COLOR = QtGui.QColor(0, 250, 20, alpha=200)  # Green
BAD_PLUGSITE_COLOR = QtGui.QColor(255, 0, 0, alpha=255)  # Red
UNPLUGED_COLOR = QtGui.QColor(0, 20, 250, alpha=200)  # Blue


class CellItem(QtGui.QGraphicsItem):
    """
    This is the parent item containing all the objects within the cell
    """

    def __init__(self, metaphase, parent=None):
        QtGui.QGraphicsItem.__init__(self, parent=parent)

        self.N = int(metaphase.KD.params['N'])
        self.Mk = int(metaphase.KD.params['Mk'])
        self.mt = metaphase  # Metaphase instance

        self.spbR = SPBItem(0, parent=self)
        self.spbL = SPBItem(-1, parent=self)

        self.cens_A = []
        self.cens_B = []

        for n in range(self.N):
            self.cens_A.append(CentromereItem(n, -1, parent=self))
            self.cens_B.append(CentromereItem(n, 0, parent=self))

        self.time_point = -1
        self.gotoTime(0)

    def advance(self):
        self.time_point += 1
        return self.gotoTime(self.time_point)

    def gotoTime(self, time_point):
        if time_point >= self.mt.KD.num_steps - 1:
            return False

        self.time_point = time_point

        self.spbR.gotoTime(self.time_point)
        self.spbL.gotoTime(self.time_point)

        for n in range(self.N):
            cen_A = self.cens_A[n]
            cen_B = self.cens_B[n]

            cen_A.gotoTime(self.time_point)
            cen_B.gotoTime(self.time_point)

    def boundingRect(self):
        """
        As this is not easily changed (or not supposed to, we give a
        large box
        """
        height = self.N * ( self.Mk * PLUGSITE_HEIGHT) * 5
        width = self.mt.KD.spbR.traj[-1] * 6
        return QtCore.QRectF(-width / 2, - height / 2, width, height)

    def paint(self, painter, option, widget):
        painter.setBrush(QtCore.Qt.white)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.1))
        painter.drawRoundedRect(self.boundingRect(), 30, 100,
            QtCore.Qt.RelativeSize)


class CentromereItem(QtGui.QGraphicsItem):

    def __init__(self, n, side, parent=None):

        QtGui.QGraphicsItem.__init__(self, parent=parent)
        self.graphcell = parent
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)

        N = int(self.graphcell.mt.KD.params['N'])
        Mk = int(self.graphcell.mt.KD.params['Mk'])

        self.n = n
        self.side = side
        self.ch = self.graphcell.mt.KD.chromosomes[n]
        if side == 0:
            self.traj = self.ch.cen_A.traj
        else:
            self.traj = self.ch.cen_B.traj

        self.y = (n - (N - 1) / 2.) * 0.55

        self.width = PLUGSITE_WIDTH * 2
        self.height = Mk * PLUGSITE_HEIGHT + (Mk -1) * PLUGSITE_OFFSET
        self.setZValue(n)

        self.plugsites = []
        for m in range(Mk):
            if side == 0:
                plugsite = self.ch.cen_A.plugsites[m]
            else:
                plugsite = self.ch.cen_B.plugsites[m]
            self.plugsites.append(PlugSiteItem(m, plugsite, self, parent))

    def gotoTime(self, time_point):
        try:
            x = self.traj[time_point]
        except IndexError:
            x = self.traj[-1]
        pos = QtCore.QPointF(x, self.y)
        self.setPos(pos)

        # Move PlugSite
        for p in self.plugsites:
            p.gotoTime(time_point)

    def shape(self):
        self.path = QtGui.QPainterPath()
        self.path.addEllipse(self.width, self.height,
                             2 * self.width, 2 * self.height)
        return self.path

    def paint(self, painter, option, widget):
        if self.ch.cen_A.is_attached() and self.ch.cen_B.is_attached():
            brush = QtGui.QBrush(CH_COLOR)
        else:
            brush = QtGui.QBrush(ACTIVE_SAC_COLOR)
        painter.setBrush(brush)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.01))
        center = self.pos()
        painter.drawEllipse(center, self.width, self.height)

    def boundingRect(self):
        adjust = 1.
        return QtCore.QRectF(-self.width - adjust, - self.height - adjust,
                             2 * self.width + adjust, 2 * self.height + adjust)


class PlugSiteItem(QtGui.QGraphicsItem):

    def __init__(self, m, plugsite, kineto, parent=None):

        QtGui.QGraphicsItem.__init__(self, parent=parent)
        self.kineto = kineto

        Mk = int(self.kineto.graphcell.mt.KD.params['Mk'])
        self.m = m
        self.sim = plugsite

        self.width = PLUGSITE_WIDTH
        self.height = PLUGSITE_HEIGHT

        self.y = self.kineto.y + (m - (Mk - 1) / 2.) * 0.11
        self.setZValue(self.kineto.n * (1 + m))

        self.gotoTime(0)

    def gotoTime(self, time_point):
        try:
            x = self.sim.traj[time_point]
        except IndexError:
            x = self.sim.traj[-1]
        pos = QtCore.QPointF(x, self.y)
        self.setPos(pos)

        # Change color according to state
        self.setState(time_point)

    def setState(self, time_point=0):
        """
        Change color according to state at specific time_point
        """
        if self.sim.state_hist[time_point] == 0:
            brush = QtGui.QBrush(UNPLUGED_COLOR)
        elif self.sim.is_correct(time_point):
            brush = QtGui.QBrush(GOOD_PLUGSITE_COLOR)
        else:
            brush = QtGui.QBrush(BAD_PLUGSITE_COLOR)
        self.color = brush

    def shape(self):
        self.path = QtGui.QPainterPath()
        self.path.addEllipse(-self.width / 2., - self.height / 2,
                             self.width, self.height)
        return self.path

    def paint(self, painter, option, widget):
        brush = self.color
        painter.setBrush(brush)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.01))
        center = self.pos()
        painter.drawEllipse(center, self.width, self.height)

    def boundingRect(self):
        adjust = 1.
        return QtCore.QRectF(-self.width - adjust, - self.height - adjust,
                             2 * self.width + adjust, 2 * self.height + adjust)


class SPBItem(QtGui.QGraphicsItem):

    def __init__(self, side, parent=None):
        QtGui.QGraphicsItem.__init__(self, parent=parent)
        self.graphcell = parent
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
        if side == 0:
            self.traj = self.graphcell.mt.KD.spbR.traj
        else:
            self.traj = self.graphcell.mt.KD.spbL.traj

        self.width = SPB_WIDTH
        self.height = SPB_HEIGHT

    def gotoTime(self, time_point):
        try:
            x = self.traj[time_point]
        except IndexError:
            x = self.traj[-1]
        y = 0.
        pos = QtCore.QPointF(x, y)
        self.setPos(pos)

    def shape(self):
        self.path = QtGui.QPainterPath()
        self.path.addEllipse(-self.width / 2., - self.height / 2,
                             self.width, self.height)
        return self.path

    def paint(self, painter, option, widget):
        brush = QtGui.QBrush(SPB_COLOR)
        painter.setBrush(brush)
        painter.setPen(QtGui.QPen(QtCore.Qt.black, 0.01))
        center = self.pos()
        painter.drawEllipse(center, self.width, self.height)

    def boundingRect(self):
        adjust = 1.
        return QtCore.QRectF(-self.width - adjust, - self.height - adjust,
                             2 * self.width + adjust, 2 * self.height + adjust)

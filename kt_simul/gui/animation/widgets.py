#-*- coding: utf-8 -*-

from PySide import QtCore, QtGui

import math
import os

from kt_simul.gui.animation.items import CellItem

FRAME_RATE = 25  # image/seconds


class InteractiveCellWidget(QtGui.QWidget):
    """
    """

    def __init__(self, metaphase):
        """
        """
        super(InteractiveCellWidget, self).__init__()

        self.timerId = 0
        self.meta = metaphase

        self.setWindowTitle("Kinetochore simulation")

        self.ccw = ControlCellWidget(metaphase)

        # This widget will contain self.view and self.textArea
        self.mainWidget = QtGui.QWidget()
        self.view = ViewCellWidget(metaphase)
        self.textArea = QtGui.QPlainTextEdit()
        self.textArea.setReadOnly(True)

        self.connect(self.ccw.playButton, QtCore.SIGNAL('clicked()'),
                     self.play)
        self.connect(self.ccw.pauseButton, QtCore.SIGNAL('clicked()'),
                     self.pause)
        self.ccw.slider.valueChanged.connect(self.gotoTime)

        # Set layouts
        hbox = QtGui.QHBoxLayout()
        hbox.setSpacing(5)
        self.mainWidget.setLayout(hbox)
        hbox.addWidget(self.view)
        hbox.addWidget(self.textArea)

        vbox = QtGui.QVBoxLayout()
        vbox.setSpacing(5)
        vbox.addWidget(self.mainWidget)
        vbox.addWidget(self.ccw)
        self.setLayout(vbox)

        # Print report for time 0
        report = self.meta.get_report()
        self.textArea.setPlainText(report)

    def pause(self):
        if not self.timerId:
            pass
        else:
            self.killTimer(self.timerId)
            self.timerId = 0

    def gotoTime(self, time):
        if self.ccw.slider.hasFocus():
            #self.pause()
            rects = [self.view.cell.boundingRect()]
            self.view.updateScene(rects)
            self.view.cell.gotoTime(time)

            # Print report
            report = self.meta.get_report(int(time))
            self.textArea.setPlainText(report)
            #self.timerId = 0

    def play(self):
        # I'll re-implement this when user can move things in graphcell
        if not self.timerId:
            self.timerId = self.startTimer(1000 / FRAME_RATE)

    def startAnim(self):
        if not self.timerId:
            self.timerId = self.startTimer(1000 / FRAME_RATE)

    def timerEvent(self, event):
        itemsMoved = False
        rects = [self.view.cell.boundingRect()]
        self.view.updateScene(rects)
        if self.view.cell.advance():
            itemsMoved = True
            self.ccw.slider.setValue(self.view.cell.time_point)
        if not itemsMoved:
            self.killTimer(self.timerId)
            self.timerId = 0


class ControlCellWidget(QtGui.QWidget):

    def __init__(self, metaphase, parent=None):
        super(ControlCellWidget, self).__init__(parent)
        numsteps = int(metaphase.KD.params['span'] / metaphase.KD.params['dt'])

        orientation = QtCore.Qt.Horizontal
        self.slider = QtGui.QSlider(orientation)
        self.slider.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider.setTickPosition(QtGui.QSlider.TicksBothSides)
        self.slider.setTickInterval(10)
        self.slider.setSingleStep(1)

        self.slider.setMinimum(0)
        self.slider.setMaximum(numsteps - 1)
        self.slider.setValue(0)

        self.playButton = QtGui.QPushButton(self)
        play_icon = os.path.join(os.path.dirname(__file__),
                                 "images", "player_play.svg")
        self.playButton.setIcon(QtGui.QIcon(play_icon))

        self.pauseButton = QtGui.QPushButton(self)
        pause_icon = os.path.join(os.path.dirname(__file__),
                                 "images", "player_pause.svg")
        self.pauseButton.setIcon(QtGui.QIcon(pause_icon))

        hbox = QtGui.QHBoxLayout()
        hbox.setSpacing(5)
        hbox.addWidget(self.playButton)
        hbox.addWidget(self.pauseButton)
        hbox.addWidget(self.slider)
        self.setLayout(hbox)


class ViewCellWidget(QtGui.QGraphicsView):

    def __init__(self, metaphase):
        super(ViewCellWidget, self).__init__()

        self.setRenderHint(QtGui.QPainter.Antialiasing)

        self.timerId = 0
        self.cell = CellItem(metaphase)
        scene = QtGui.QGraphicsScene(self)

        self.setCacheMode(QtGui.QGraphicsView.CacheBackground)
        #self.setViewportUpdateMode(QtGui.QGraphicsView.BoundingRectViewportUpdate)
        self.setRenderHint(QtGui.QPainter.Antialiasing)
        self.setTransformationAnchor(QtGui.QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QtGui.QGraphicsView.AnchorViewCenter)

        scene.setSceneRect(-9, -4, 18, 8)
        self.setScene(scene)
        scene.addItem(self.cell)
        self.scale(40, 40)

    def wheelEvent(self, event):
        self.scaleView(math.pow(2.0, event.delta() / 240.0))

    def scaleView(self, scaleFactor):
        factor = self.matrix().scale(scaleFactor, scaleFactor).mapRect(
            QtCore.QRectF(0, 0, 1, 1)).width()
        if factor < 0.1 or factor > 200:
            return
        self.scale(scaleFactor, scaleFactor)

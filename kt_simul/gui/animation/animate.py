"""
Run and play simulation animation in 2D (need PySide)
"""

import sys
import logging

from PyQt4 import QtGui, QtCore

from kt_simul.gui.animation import InteractiveCellWidget

log = logging.getLogger(__name__)


class Animator:
    """
    """

    def __init__(self, meta_instance):
        """
        """
        self.meta = meta_instance

    def play(self,):
        """
        """

        log.info("Playing animation")

        self._init_qt()
        QtCore.qsrand(QtCore.QTime(0, 0, 0).secsTo(QtCore.QTime.currentTime()))
        widget = InteractiveCellWidget(self.meta)
        self._create_window(widget)

    def _init_qt(self):
        """
        """
        self.app_created = False
        self.app = QtCore.QCoreApplication.instance()
        if self.app is None:
            self.app = QtGui.QApplication(sys.argv)
            self.app_created = True
        self.app.references = set()

    def _create_window(self, window):
        """
        Create a QT window in Python, or interactively in IPython with QT GUI
        event loop integration.
        """
        self.app.references.add(window)
        window.show()
        if self.app_created:
            self.app.exec_()

# -*- coding: utf-8 -*-

"""
Graphical User Interface for the kinetochore dynamics simulation
"""

import sys
import os
import time

from numpy import mod

from PySide import QtCore, QtGui

import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'
from kt_simul.io.xml_handler import ParamTree

from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

from kt_simul.core.simul_spindle import MEASUREFILE, PARAMFILE

from kt_simul.gui.param_seter import SetParameters, SetMeasures
from kt_simul.gui.game import InteractiveCellWidget
from kt_simul.gui.canvas import MyMplCanvas
from kt_simul.gui.sig_metaphase import SigMetaphase

__all__ = ['MainWindow']


class MainWindow(QtGui.QMainWindow):

    def __init__(self, parent=None):

        self.paramtree = ParamTree(PARAMFILE)
        self.measuretree = ParamTree(MEASUREFILE, adimentionalized=False)
        self.measures = self.measuretree.absolute_dic

        self.mt = None

        QtGui.QMainWindow.__init__(self, parent)
        self.centralwidget = QtGui.QWidget()

        self.setCentralWidget(self.centralwidget)
        self.create_docks()
        self.create_tabs()
        self.create_buttons()

        #self.prepare_simulation()

    def create_docks(self):

        #Parameter Setting in a Dock Widget
        self.paramdock = QtGui.QDockWidget('Parameters Setting')
        self.setParameters = SetParameters(self.paramtree)
        scrollArea = QtGui.QScrollArea()
        scrollArea.setWidget(self.setParameters)
        self.paramdock.setWidget(scrollArea)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.paramdock)

        #Measures Setting in another Dock Widget
        self.measuredock = QtGui.QDockWidget('Measures Setting')
        self.setMeasures = SetMeasures(self.measuretree)
        scrollArea = QtGui.QScrollArea()
        scrollArea.setWidget(self.setMeasures)
        self.measuredock.setWidget(scrollArea)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.measuredock)

    def create_tabs(self):
        #Text Area
        s = ("Welcome to the S.Pombe kinetochore motion "
             "Simulator")
        self.simLog = QtGui.QTextEdit(s)
        self.simLog.setReadOnly(True)

        #Plotting Areas
        span = self.paramtree.relative_dic['span']

        self.plotarea1 = MyMplCanvas(span)
        mpl_toolbar = NavigationToolbar(self.plotarea1, self.centralwidget)
        vbox1 = QtGui.QVBoxLayout()
        vbox1.setSpacing(5)
        vbox1.addWidget(self.plotarea1)
        vbox1.addWidget(mpl_toolbar)
        self.w1 = QtGui.QWidget()
        self.w1.setLayout(vbox1)

        self.plotarea2 = MyMplCanvas(span)
        mpl_toolbar = NavigationToolbar(self.plotarea2, self.centralwidget)
        vbox2 = QtGui.QVBoxLayout()
        vbox2.setSpacing(5)
        vbox2.addWidget(self.plotarea2)
        vbox2.addWidget(mpl_toolbar)
        self.w2 = QtGui.QWidget()
        self.w2.setLayout(vbox2)

        #All this goes in a tab widget
        self.tabWidget = QtGui.QTabWidget()
        self.tabWidget.setTabsClosable(True)
        self.tabWidget.addTab(self.simLog, "Log")
        #self.tabWidget.tabCloseRequested.connect(self.tabWidget.removeTab)

        self.connect(self.tabWidget, QtCore.SIGNAL('tabCloseRequested(int)'),
                     self.closeTab)
        #              self.closeTab)
        #self.tabWidget.addTab(w1, "All trajectories")
        #self.tabWidget.addTab(w2, "One trajectory")

    def closeTab(self, idx):
        self.tabWidget.removeTab(idx)

    def create_buttons(self):
        #Buttons
        #self.buttonGroup = QtGui.QButtonGroup()
        runButton = QtGui.QPushButton('Run the simulation')
        exec_icon = os.path.join(os.path.dirname(__file__),
                                 "images", "exec.svg")
        runButton.setIcon(QtGui.QIcon(exec_icon))
        self.connect(runButton, QtCore.SIGNAL('clicked()'), self.run_simulation)

        showTrajButton = QtGui.QPushButton('Show trajectories')
        self.connect(showTrajButton, QtCore.SIGNAL('clicked()'), self.show_trajs)

        self.traj_num = 0
        showOneButton = QtGui.QPushButton('Show one trajectory')
        self.connect(showOneButton, QtCore.SIGNAL('clicked()'), self.show_one)

        self.interactiveButton = QtGui.QRadioButton('Interactive Simulation')
        self.interactiveButton.setChecked(True)
        #self.buttonGroup.addButton(runButton)
        #Progress Bar
        self.progressBar = QtGui.QProgressBar()

        hbox = QtGui.QHBoxLayout()
        hbox.setSpacing(5)
        hbox.addWidget(runButton)
        hbox.addWidget(showTrajButton)
        hbox.addWidget(showOneButton)
        hbox.addWidget(self.progressBar)
        hbox.addWidget(self.interactiveButton)

        vbox = QtGui.QVBoxLayout()
        vbox.setSpacing(5)
        vbox.addLayout(hbox)  #self.buttonGroup)
        vbox.addWidget(self.tabWidget)

        self.centralwidget.setLayout(vbox)

        self.createActions()
        self.createMenus()

        self.createToolBars()
        self.createStatusBar()
        currentDirName = os.path.abspath(os.path.curdir)

        currentFileName = '_'.join(time.asctime().split()[:-1])
        currentFileName = 'simul_' + \
                            '_'.join(currentFileName.split(':')) + \
                            '.xml'
        currentFileName = os.path.join(currentDirName, currentFileName)

        self.setCurrentFile(currentFileName)
        self.setWindowTitle(self.tr("Kinetochore Dynamics Simulation"))
        self.setMinimumSize(160, 160)
        self.resize(1000, 600)

    def prepare_simulation(self):

        plug_idx = self.attachCombo.currentIndex()
        initial_plug = self.attachment_list[plug_idx]

        self.mt = SigMetaphase(self.paramtree, self.measuretree,
                               initial_plug=initial_plug)

        self.progressBar.setMaximum(int(self.paramtree.absolute_dic['span']))

        self.progressBar.setMinimum(0)
        self.connect(self.mt, QtCore.SIGNAL('plugCheckPoint'),
                     self.active_checkpoint)
        self.connect(self.mt,  QtCore.SIGNAL('stepDone'),
                     self.progressBar.setValue)
        self.connect(self.mt,  QtCore.SIGNAL('simulDone'),
                     self.print_report)

        run_interactively = self.interactiveButton.isChecked()

        if run_interactively:
            self.iw = InteractiveCellWidget(self.mt)

            idx = self.tabWidget.currentIndex() + 1
            self.tabWidget.insertTab(idx, self.iw, "Interactive Simulation")
            self.tabWidget.setCurrentIndex(idx)

            self.connect(self.mt,  QtCore.SIGNAL('simulDone'),
                         self.iw.startAnim)

    def run_simulation(self):

        self.prepare_simulation()
        QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        self.mt.sig_simul()
        QtGui.QApplication.restoreOverrideCursor()

    def show_trajs(self):

        self.plotarea1.update_figure(self.mt)
        self.tabWidget.insertTab(1, self.w1, "All trajectories")
        self.tabWidget.setCurrentIndex(1)

    def show_one(self):

        #self.tabWidget.removeTab(2)
        t_num = mod(self.traj_num, 3)
        self.plotarea2.update_figure(self.mt, t_num)
        self.traj_num += 1
        self.tabWidget.insertTab(2, self.w2, "Trajectory of chromosome # %i" %(t_num + 1) )
        self.tabWidget.setCurrentIndex(2)

    def update_progressBar(self, val):
        self.progressBar.setValue(val)

    def print_report(self, report):
        self.simLog.append("Simulation's done!")
        for l in report:
            ls = l
            self.simLog.append(ls)

        self.progressBar.setValue(0)

    def active_checkpoint(self, cp):
        if cp:
            self.statusBar().showMessage(self.tr("Active Plug/unplug checkpoint"), 2000)
        else:
            self.statusBar().showMessage(self.tr(" "), 2000)

    def closeEvent(self, event):
        if self.maybeSave():
            event.accept()
        else:
            event.ignore()

    def newFile(self):
        if self.maybeSave():
            self.textEdit.clear()
            self.setCurrentFile("")

    def open(self):
        if self.maybeSave():
            fileName, selectedFilter = QtGui.QFileDialog.getOpenFileName(self, filter="XML files (*.xml);;All Files (*.*)")
            self.loadFile(fileName)

        self.setCurrentFile(fileName)

    def save(self):
        # if not os.path.isfile(self.curFile):
        #     self.curFile = default_filename
        #     #return self.saveAs()
        # else:
        return self.saveFile(self.curFile)

    # TODO: The saveAs functions trigs a segfault
    # def saveAs(self):
    #     #'test_save'#
    #     saveDiag = QtGui.QFileDialog()
    #     fileName = saveDiag.getSaveFileName(self,
    #                                         self.tr("Save Simulation"),
    #                                         self.tr("default.xml"),
    #                                         self.tr("Simulations (*.xml)"))
    ### SegFaults here, I don't know why ###
    #     print filename
    #     if fileName.isEmpty():
    #         fileName = "default.xml"

    #     self.saveFile(fileName)
    #     del saveDiag

    def about(self):
        QtGui.QMessageBox.about(self, self.tr("About Application"),
            self.tr("This little stuff allows to simulate "
                    "kinetochore dynamics in S.Pombe "))

    def documentWasModified(self):
        self.setWindowModified(self.setParameters.isModified)

    def createActions(self):

        open_icon = os.path.join(os.path.dirname(__file__),
                                 "images", "open.png")
        self.openAct = QtGui.QAction(QtGui.QIcon(open_icon),
                                     self.tr("&Open Simulation File"), self)
        self.openAct.setShortcut(self.tr("Ctrl+O"))
        self.openAct.setStatusTip(self.tr("Open an existing simulation file"))
        self.connect(self.openAct, QtCore.SIGNAL("triggered()"), self.open)
        save_icon = os.path.join(os.path.dirname(__file__),
                                 "images", "save.png")
        self.saveAct = QtGui.QAction(QtGui.QIcon(save_icon),
                                     self.tr("&Save"), self)
        self.saveAct.setShortcut(self.tr("Ctrl+S"))
        self.saveAct.setStatusTip(self.tr("Save the simulation to disk"))
        self.connect(self.saveAct, QtCore.SIGNAL("triggered()"), self.save)

        # self.saveAsAct = QtGui.QAction(self.tr("Save &As..."), self)
        # self.saveAsAct.setStatusTip(self.tr("Save the simulation under a new name"))
        # self.connect(self.saveAsAct, QtCore.SIGNAL("triggered()"), self.saveAs)

        self.exitAct = QtGui.QAction(self.tr("E&xit"), self)
        self.exitAct.setShortcut(self.tr("Ctrl+Q"))
        self.exitAct.setStatusTip(self.tr("Exit the application"))
        self.connect(self.exitAct, QtCore.SIGNAL("triggered()"), self, QtCore.SLOT("close()"))

        self.aboutAct = QtGui.QAction(self.tr("&About"), self)
        self.aboutAct.setStatusTip(self.tr("Show the application's About box"))
        self.connect(self.aboutAct, QtCore.SIGNAL("triggered()"), self.about)

        self.aboutQtAct = QtGui.QAction(self.tr("About &Qt"), self)
        self.aboutQtAct.setStatusTip(self.tr("Show the Qt library's About box"))
        self.connect(self.aboutQtAct, QtCore.SIGNAL("triggered()"), QtGui.qApp, QtCore.SLOT("aboutQt()"))

    def createMenus(self):
        self.fileMenu = self.menuBar().addMenu(self.tr("&File"))
        self.fileMenu.addAction(self.openAct)
        self.fileMenu.addAction(self.saveAct)
        # self.fileMenu.addAction(self.saveAsAct)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAct)

        self.menuBar().addSeparator()

        self.helpMenu = self.menuBar().addMenu(self.tr("&Help"))
        self.helpMenu.addAction(self.aboutAct)
        self.helpMenu.addAction(self.aboutQtAct)

    def createToolBars(self):
        self.fileToolBar = self.addToolBar(self.tr("File"))
        self.fileToolBar.addAction(self.openAct)
        self.fileToolBar.addAction(self.saveAct)

        self.attachCombo = QtGui.QComboBox()
        self.attachCombo.addItem("No attachment")  #0
        self.attachCombo.addItem("Amphitelic attachment") #1
        self.attachCombo.addItem("Merotelic attachment")  #2
        self.attachCombo.addItem("Syntelic attachment")  #3
        self.attachCombo.addItem("Monotelic attachment")  #4
        self.attachCombo.addItem("Random attachment")  #5

        self.configToolBar = self.addToolBar(self.tr("config"))
        self.configToolBar.addWidget(self.attachCombo)

        self.attachment_list = ['null',
                                'amphitelic',
                                'merotelic',
                                'syntelic',
                                'monotelic',
                                'random']
        self.attachCombo.setCurrentIndex(5)

    def createStatusBar(self):
        self.statusBar().showMessage(self.tr("Ready"))

    def maybeSave(self):
        if self.setParameters.isModified:
            ret = QtGui.QMessageBox.warning(self, self.tr("Application"),
                        self.tr("The document has been modified.\n"
                                "Do you want to save your changes?"),
                        QtGui.QMessageBox.Yes | QtGui.QMessageBox.Default,
                        QtGui.QMessageBox.No,
                        QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Escape)
            if ret == QtGui.QMessageBox.Yes:
                return self.save()
            elif ret == QtGui.QMessageBox.Cancel:
                return False
        return True

    def loadFile(self, fileName):

        QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        self.prepare_simulation()
        tmp_m = get_fromfile(str(fileName))
        self.paramtree = tmp_m.paramtree
        self.measuretree = tmp_m.measuretree

        self.removeDockWidget(self.paramdock)
        del self.paramdock
        self.removeDockWidget(self.measuredock)
        del self.measuredock
        self.create_docks()
        self.prepare_simulation()
        self.mt = SigMetaphase(tmp_m.paramtree, tmp_m.measuretree)
        self.mt.KD = tmp_m.KD
        del tmp_m
        self.mt.emit(QtCore.SIGNAL('simulDone'),self.mt.report)

        QtGui.QApplication.restoreOverrideCursor()

        self.setCurrentFile(fileName)
        self.statusBar().showMessage(self.tr("Loaded"), 2000)


    def saveFile(self, fileName):

        if not fileName.endsWith('.xml'):
            xmlfname = fileName+'.xml'
            datafname = fileName+'_data.npy'
        else:
            xmlfname = fileName
            datafname = fileName.split('.')[-2]+'_data.npy'

        self.mt.write_results(str(xmlfname), str(datafname))

        self.setCurrentFile(fileName);
        self.statusBar().showMessage(self.tr("File saved"), 2000)
        return True

    def setCurrentFile(self, fileName):

        self.curFile = fileName
        self.setParameters.setModified(False)
        self.setWindowModified(False)

    def strippedName(self, fullFileName):
        return QtCore.QFileInfo(fullFileName).fileName()

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)

    mainwindow = MainWindow()
    mainwindow.show()
    sys.exit(app.exec_())

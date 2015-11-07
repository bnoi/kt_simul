import numpy as np

import vispy
from vispy import geometry
from vispy import scene

from . import StructureViewer
from . import StructureWidget

__class__ = ["SpindleWidget", "SpindleViewer"]


class SpindleViewer(StructureViewer):
    """This viewer only use Vispy and automatic backend associated.
    """

    def __init__(self, metaphase=None):

        super().__init__()

        if metaphase:
            self.add_structure(metaphase)

    def move(self, time_point):

        if time_point == self.duration:
            time_point = 0

        dt = self.metaphase.params['dt']
        mess = "Time : {:.0f}/{:.0f} ({:.0f}s/{:.0f}s)".format(time_point + 1, self.duration + 1,
                                                               time_point * dt, self.duration * dt)
        self.set_status(mess, 'Time')

        def _move(sphere, point):
            sphere.transform.reset()
            sphere.transform.translate(point.traj.loc[time_point, self.coords].values)

        # Move poles
        _move(*self.points[self.spbL_idx])
        _move(*self.points[self.spbR_idx])

        # Move centromeres
        for cen_A_idx, cen_B_idx in self.ch_idxs:
            _move(*self.points[cen_A_idx])
            _move(*self.points[cen_B_idx])

        self.time_point = time_point

    def add_structure(self, metaphase):

        self.metaphase = metaphase
        self.structure = self.metaphase.spindle
        self.duration = self.structure.point_hist.keys()[-1]

        # Get points indexes
        self.spbL_idx = self.metaphase.spindle.spbL.idx
        self.spbR_idx = self.metaphase.spindle.spbR.idx
        self.ch_idxs = []
        for ch in self.structure.chromosomes:
            self.ch_idxs.append((ch.cen_A.idx, ch.cen_B.idx))

        # Create spheres
        self.structure.points[self.spbL_idx].color = "gray"
        self.structure.points[self.spbR_idx].color = "gray"
        self.create_sphere(self.structure.points[self.spbL_idx], time_point=0, radius=0.8)
        self.create_sphere(self.structure.points[self.spbR_idx], time_point=0, radius=0.8)

        cmap = vispy.color.get_colormap('husl', value=0.5)
        colors = cmap.map(np.linspace(0, 1, len(self.ch_idxs)))
        for color, (cen_A_idx, cen_B_idx) in zip(colors, self.ch_idxs):

            self.structure.points[cen_A_idx].color = color
            self.structure.points[cen_B_idx].color = color

            self.create_sphere(self.structure.points[cen_A_idx], time_point=0, radius=0.8)
            self.create_sphere(self.structure.points[cen_B_idx], time_point=0, radius=0.8)

        # Set a confortable zoom level
        # TODO: automatic scale factor calculation
        self.view.camera.scale_factor = 10

        # Create border cell
        self.create_border_cell()

        # Setup time
        self.move(time_point=0)

    def create_border_cell(self):

        mdata = geometry.create_sphere(64, 64,  radius=15)
        borders = scene.visuals.Mesh(meshdata=mdata, shading='flat', color="#33333333")
        self.view.add(borders)


class SpindleWidget(StructureWidget):
    """This viewer use Vispy and Qt to display some controls and informations about the simulation
    """

    def __init__(self):

        StructureWidget.QWidget.__init__(self)

        self.setMinimumSize(800, 600)
        self.setWindowTitle('Kymo like')

        # Config time slider
        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.slider.valueChanged.connect(self.move)
        self.slider.setEnabled(False)

        # Config play button
        self.play_button = QtGui.QPushButton('Start', self)
        self.play_button.setEnabled(False)
        self.play_button.clicked.connect(self.play_simu)

        # Add multiples status bar
        status_title = ["FPS", "Time"]
        self.status = {title: QtGui.QStatusBar(parent=self) for title in status_title}

        # Status layout
        self.status_widget = QtGui.QWidget()
        self.status_layout = QtGui.QHBoxLayout()
        self.status_widget.setLayout(self.status_layout)
        for i, title in enumerate(status_title):
            self.status_layout.addWidget(self.status[title], i+1)
            self.status[title].showMessage("{} : ".format(title))
        self.status_layout.addItem(QtGui.QSpacerItem(100, 0))

        # Play layout
        self.play_widget = QtGui.QWidget()
        self.play_layout = QtGui.QHBoxLayout()
        self.play_widget.setLayout(self.play_layout)

        self.play_layout.addWidget(self.play_button, 0)
        self.play_layout.addWidget(self.slider, 1)

        # Screen layout
        self.screen_widget = QtGui.QWidget()
        self.screen_layout = QtGui.QHBoxLayout()
        self.screen_widget.setLayout(self.screen_layout)

        self.screen_layout.addWidget(self.canvas.native, 0)
        # self.screen_layout.addWidget(self.slider, 1)

        # Main Layout
        self.layout = QtGui.QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.addWidget(self.screen_widget, 0)
        self.layout.addItem(QtGui.QSpacerItem(100, 0))
        self.layout.addWidget(self.play_widget, 2)
        self.layout.addWidget(self.status_widget, 3)

        # Setup timer
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.update_simu)
        self.simu_playing = False

    def set_status(self, message, status_type=None):
        self.status[status_type].showMessage(message)

    def show_fps(self, fps):
        mess = "FPS : {:.1f}".format(fps)
        self.set_status(mess, 'FPS')

    def add_spindle(self, metaphase):

        SimuViewer.add_spindle(self, metaphase)

        # Set time slider and others widgets
        self.slider.setRange(0, self.spindle.duration)
        self.slider.setSingleStep(self.spindle.dt)
        self.slider.setEnabled(True)

        self.play_button.setEnabled(True)

    def play_simu(self):
        """
        """
        if not self.simu_playing:
            self.play_button.setText("Pause")
            self.timer.start(100)
            self.simu_playing = True
        else:
            self.play_button.setText("Start")
            self.timer.stop()
            self.simu_playing = False

    def update_simu(self):
        self.slider.setValue(self.time_point + 1)

    def show(self):

        QtGui.QWidget.show(self)

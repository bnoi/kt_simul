from collections import OrderedDict

from vispy import geometry
from vispy import scene
from vispy import app
from vispy import visuals

qtapp = app.use_app('pyqt4')

QtCore = qtapp.backend_module.QtCore
QtGui = qtapp.backend_module.QtGui

__class__ = ["StructureWidget", "StructureViewer"]


class StructureCanvas(scene.SceneCanvas):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def on_close(self, event):
        app.quit()
        super().on_close(event)


class StructureViewer():
    """This viewer only use Vispy and automatic backend associated.

    Usage
    -----

    >>> v = StructureViewer()
    >>> v.add_structure(sprg)
    >>> v.play()
    """

    def __init__(self, structure=None):

        self.canvas = StructureCanvas(keys='interactive', resizable=True,
                                      fullscreen=False, bgcolor="#fffaf0", show=True)

        self.canvas.measure_fps(0.1, self.show_fps)

        self.view = self.canvas.central_widget.add_view()
        self.view.camera = scene.cameras.TurntableCamera()

        self.points = {}
        self.time_point = 0
        self.coords = ['x', 'y', 'z']
        self.save = None
        self._has_been_saved = False

        self.create_axes()
        self.switch_axes()

        # Setup status
        self.status = scene.visuals.Text('', anchor_x='left', anchor_y='top',
                                         color='black', parent=self.view)
        self.status.pos = self.canvas.size[0] - self.canvas.size[0] * 0.98, \
                          self.canvas.size[1] - self.canvas.size[1] * 0.98

        self.status_messages = OrderedDict([('FPS', ''), ('Time', '')])

        if structure:
            self.add_structure(structure)

    def set_status(self, message, status_type):

        self.status_messages[status_type] = message
        self.status.text = ' | '.join(['{}'.format(v) for k, v in self.status_messages.items()])

    def show_fps(self, fps):

        mess = "FPS : {:.1f}".format(fps)
        self.set_status(mess, 'FPS')

    def create_axes(self):

        self.axes = {}
        coords = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        names = ["x", "y", "z"]
        colors = ["red", "green", "blue"]

        for coord, color, name in zip(coords, colors, names):

            mdata = geometry.create_cylinder(128, 128, radius=[0.1, 0.1],
                                             length=5, offset=False)
            axis = scene.visuals.Mesh(meshdata=mdata, shading='flat', color=color)
            t = visuals.transforms.MatrixTransform()
            t.rotate(90, coord)
            axis.transform = t

            self.view.add(axis)
            self.axes[name] = axis
            axis.visible = False

    def switch_axes(self):

        for name, axis in self.axes.items():
            axis.visible = not axis.visible

    def create_sphere(self, point, time_point, radius):

        mdata = geometry.create_sphere(64, 64,  radius=radius)

        color = point.color
        if isinstance(color, type(None)):
            color = "gray"

        sphere = scene.visuals.Mesh(meshdata=mdata, shading='flat', color=color)

        t = visuals.transforms.MatrixTransform()
        sphere.transform = t
        sphere.transform.translate(point.traj.loc[time_point, self.coords].values)

        self.points[point.idx] = (sphere, point)
        self.view.add(sphere)
        return sphere

    def move(self, time_point):

        if time_point == self.duration + 1:
            time_point = 0

        mess = "Time point : {:.0f}/{:.0f}".format(time_point, self.duration)
        self.set_status(mess, 'Time')

        for idx, (sphere, point) in self.points.items():
            sphere.transform.reset()
            sphere.transform.translate(point.traj.loc[time_point, self.coords].values)

        self.time_point = time_point

    def add_structure(self, structure):

        self.structure = structure
        self.duration = self.structure.point_hist.keys()[-1]

        for idx, point in self.structure.points.items():
            self.create_sphere(point, time_point=0, radius=0.8)

        # Set a confortable zoom level
        # TODO: automatic scale factor calculation
        self.view.camera.scale_factor = 10

        # Setup time
        self.move(time_point=0)

    def show(self):

        app.run()

    def quit(self):

        app.quit()

    def play(self):

        self.move(time_point=0)
        self.timer = app.Timer(connect=self.update)
        self.timer.start(0.1)

    def update(self, event):

        self.move(self.time_point + 1)


class StructureWidget(StructureViewer, QtGui.QWidget):
    """This viewer use Vispy and Qt to display some controls and informations about the simulation.

    Usage
    -----

    >>> w = StructureWidget()
    >>> w.add_structure(sprg)
    >>> w.show()
    """

    def __init__(self, structure=None):

        QtGui.QWidget.__init__(self)
        StructureViewer.__init__(self, structure)

        self.setMinimumSize(800, 600)
        self.setWindowTitle("Structure's points trajectories")

        # Config time slider
        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.slider.valueChanged.connect(self.move)
        self.slider.setEnabled(False)

        # Config play button
        self.play_button = QtGui.QPushButton('Start', self)
        self.play_button.setEnabled(False)
        self.play_button.clicked.connect(self.play)

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
        self.timer.timeout.connect(self.update)
        self.simu_playing = False

    def set_status(self, message, status_type=None):
        self.status[status_type].showMessage(message)

    def show_fps(self, fps):
        mess = "FPS : {:.1f}".format(fps)
        self.set_status(mess, 'FPS')

    def add_structure(self, structure):

        StructureViewer.add_structure(self, structure)

        # Set time slider and others widgets
        self.slider.setRange(0, self.duration)
        self.slider.setSingleStep(1)
        self.slider.setEnabled(True)

        self.play_button.setEnabled(True)

    def play(self):
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

    def update(self):
        self.slider.setValue(self.time_point + 1)

    def show(self):

        QtGui.QWidget.show(self)

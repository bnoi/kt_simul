from . import StructureViewer
from . import StructureWidget

__class__ = ["SpindleWidget", "SpindleViewer"]


class SpindleViewer(StructureViewer):
    """This viewer only use Vispy and automatic backend associated.
    """

    def __init__(self):

        StructureWidget.__init__(self)

        self.canvas = scene.SceneCanvas(keys='interactive', resizable=True,
                                        fullscreen=False, bgcolor="#fffaf0", show=True)

        self.canvas.measure_fps(0.1, self.show_fps)

        self.view = self.canvas.central_widget.add_view()
        self.view.camera = scene.cameras.TurntableCamera()

        self.points = []
        self.time_point = 0
        self.coords = ['x', 'y', 'z']

        self.create_axes()
        self.switch_axes()

        # Setup status
        self.status = scene.visuals.Text('', anchor_x='left', anchor_y='top',
                                         color='black', parent=self.view)
        self.status.pos = self.canvas.size[0] - self.canvas.size[0] * 0.98, \
                          self.canvas.size[1] - self.canvas.size[1] * 0.98

        self.status_messages = OrderedDict([('FPS', ''), ('Time', '')])

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

    def create_points(self, position, color, radius):

        mdata = geometry.create_sphere(64, 64,  radius=radius)
        point = scene.visuals.Mesh(meshdata=mdata, shading='flat', color=color)

        t = visuals.transforms.MatrixTransform()
        point.transform = t
        point.transform.translate(position)

        self.points.append(point)
        self.view.add(point)
        return point

    def move(self, time_point):

        mess = "Time : {:.0f}s/{:.0f}s".format(time_point * self.spindle.dt, self.max_t)
        self.set_status(mess, 'Time')

        # Move poles
        self.points[0].transform.reset()
        self.points[0].transform.translate(self.meta.spindle.spbL.traj.loc[time_point, self.coords].values)

        self.points[1].transform.reset()
        self.points[1].transform.translate(self.meta.spindle.spbR.traj.loc[time_point, self.coords].values)

        # Move kinetochores
        for i, ch in enumerate(self.spindle.chromosomes):

            self.points[i*2+2].transform.reset()
            self.points[i*2+2].transform.translate(ch.cen_A.traj.loc[time_point, self.coords].values)

            self.points[i*2+3].transform.reset()
            self.points[i*2+3].transform.translate(ch.cen_B.traj.loc[time_point, self.coords].values)

        self.time_point = time_point

    def add_spindle(self, metaphase):

        self.meta = metaphase
        self.spindle = self.meta.spindle
        self.max_t = self.spindle.duration / self.spindle.dt

        # Setup poles
        self.create_points(self.meta.spindle.spbL.traj.loc[0, self.coords].values, "gray", 0.8)
        self.create_points(self.meta.spindle.spbR.traj.loc[0, self.coords].values, "gray", 0.8)

        # Setup kinetochores
        cmap = vispy.color.get_colormap('husl', value=0.5)
        self.colors = cmap.map(np.linspace(0, 1, 3))
        for ch, color in zip(self.spindle.chromosomes, self.colors):
            self.create_points(ch.cen_A.traj.loc[0, self.coords].values, color, 0.5)
            self.create_points(ch.cen_B.traj.loc[0, self.coords].values, color, 0.5)

        # Create border cell
        self.create_border_cell()

        # Setup time
        self.move(time_point=0)

    def create_border_cell(self):

        mdata = geometry.create_sphere(64, 64,  radius=15)
        borders = scene.visuals.Mesh(meshdata=mdata, shading='flat', color="#33333333")

        self.view.add(borders)

    def show(self):

        app.run()

    def play_simu(self):

        self.timer = app.Timer(connect=self.update_simu)
        self.timer.start(0.1)

    def update_simu(self, event):

        self.move(self.time_point + 1)


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

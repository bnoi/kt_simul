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

        if time_point == self.duration + 1:
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

        radius = 0.5

        # Create spheres
        self.structure.points[self.spbL_idx].color = "gray"
        self.structure.points[self.spbR_idx].color = "gray"
        self.create_sphere(self.structure.points[self.spbL_idx], time_point=0, radius=radius)
        self.create_sphere(self.structure.points[self.spbR_idx], time_point=0, radius=radius)

        cmap = vispy.color.get_colormap('husl', value=0.5)
        colors = cmap.map(np.linspace(0, 1, len(self.ch_idxs)))
        for color, (cen_A_idx, cen_B_idx) in zip(colors, self.ch_idxs):

            self.structure.points[cen_A_idx].color = color
            self.structure.points[cen_B_idx].color = color

            self.create_sphere(self.structure.points[cen_A_idx], time_point=0, radius=radius)
            self.create_sphere(self.structure.points[cen_B_idx], time_point=0, radius=radius)

        # Set a confortable zoom level
        # TODO: automatic scale factor calculation
        self.view.camera.scale_factor = 10

        # Create border cell
        self.create_border_cell()

        # Setup time
        self.move(time_point=0)

    def create_border_cell(self):

        mdata = geometry.create_sphere(64, 64,  radius=15)
        borders = scene.visuals.Mesh(meshdata=mdata, shading='flat', color="#f5f5f527")
        self.view.add(borders)


class SpindleWidget(StructureWidget, SpindleViewer):
    """This viewer use Vispy and Qt to display some controls and informations about the simulation
    """

    def __init__(self, metaphase=None):

        StructureWidget.__init__(self)

        if metaphase:
            self.add_structure(metaphase)

    def add_structure(self, structure):

        SpindleViewer.add_structure(self, structure)

        # Set time slider and others widgets
        self.slider.setRange(0, self.duration)
        self.slider.setSingleStep(1)
        self.slider.setEnabled(True)

        self.play_button.setEnabled(True)

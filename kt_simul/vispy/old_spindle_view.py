"""
Simple demonstration of Box visual.
"""

from vispy import app, gloo, visuals, scene
from vispy.color import Color
from vispy.visuals.transforms import MatrixTransform


def setup(spindle):

    spindle.point_df['radius'] = 0.1
    spindle.point_df['color'] = '#a1a1a1'


    spindle.point_df.loc[spindle.spbR.idx, 'radius'] = 0.1
    spindle.point_df.loc[spindle.spbR.idx, 'color'] = '#000000'
    spindle.point_df.loc[spindle.spbL.idx, 'radius'] = 0.1
    spindle.point_df.loc[spindle.spbL.idx, 'color'] = '#000000'



    ch_colors = ['#00ffaa', '#11ffbb', '#22ffcc']
    for ch, color in zip(spindle.chromosomes, ch_colors):
        spindle.point_df.loc[ch.cen_A.idx, 'radius'] = 0.2
        spindle.point_df.loc[ch.cen_A.idx, 'color'] = color
        spindle.point_df.loc[ch.cen_B.idx, 'radius'] = 0.2
        spindle.point_df.loc[ch.cen_B.idx, 'color'] = color


class VispyOrganite():

    def __init__(self, point, spindle):
        self.point = point
        radius, color = spindle.point_df.loc[point.idx, ['radius', 'color']]
        self.sphere = scene.visuals.Sphere(radius=radius, method='ico',
                                           subdivisions=4, color=color, parent=view.scene)

        self.transform = MatrixTransform()
        self.sphere.transform = self.transform
        self.traj = point.traj * 0.1


class VispyLink():

    def __init__(self, link, spindle):
        visuals.LineVisual()

def make_anim(spindle):

    setup(spindle)

    canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)

    # Set up a viewbox to display the cube with interactive arcball
    view = canvas.central_widget.add_view()
    view.bgcolor = '#efefef'
    view.camera = 'turntable'
    view.padding = 100

    color = Color("#3f51b5")
    objects = []
    for point in spindle.points.values():
        objects.append(VispyOrganite(point, spindle))
    t = 0

    def move(app):
        global t
        if t >= spindle.point_hist.shape[0]:
            t = 0
        for obj in objects:
            obj.sphere.transform.reset()
            x, y, z = obj.traj.loc[t]
            obj.sphere.transform.translate([x, y, z])#, z))
        img = canvas.render()
        vp_io.write_png('movie/spindle_{0:03d}.png'.format(t), img)
        t += 1

    timer = app.Timer(connect=move)
    timer.start(0.1)

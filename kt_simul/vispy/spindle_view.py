"""
Simple demonstration of Box visual.
"""

from vispy import app, gloo, visuals
from vispy.geometry import create_box
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
        self.sphere = visuals.SphereVisual(radius=radius, method='ico',
                                           subdivisions=4, color=color)

        self.transform = MatrixTransform()
        self.sphere.transform = self.transform


class Canvas(app.Canvas):

    def __init__(self, spindle):
        self.sphere = visuals.SphereVisual(radius=1, method='ico',
                                           subdivisions=4, color=color)
        app.Canvas.__init__(self, keys='interactive', size=(800, 800))
        self.spindle = spindle
        self.objects = []
        for point in spindle.points.values():
            self.objects.append(VispyOrganite(point, spindle))

        self.show()
        self.t = 0

        self.timer = app.Timer(connect=self.move)
        self.timer.start(0.016)

    def move(self, event):
        if self.t >= spindle.point_hist.shape[0]:
            self.t = 0
        for obj in self.objects:
            x, y, z = self.spindle.point_hist.loc[self.t, obj.point.idx, ['x', 'y', 'z']]*10 + [400, 400, 0]
            obj.sphere.transform.reset()
            obj.sphere.transform.scale((50, 50, 0.001))
            obj.sphere.transform.translate((x, y))#, z))
        self.t += 1
        self.update()

    def on_resize(self, event):
        # Set canvas viewport and reconfigure visual transforms to match.
        vp = (0, 0, self.physical_size[0], self.physical_size[1])
        self.context.set_viewport(*vp)
        for obj in self.objects:
            obj.sphere.transforms.configure(canvas=self, viewport=vp)

    def on_draw(self, ev):
        gloo.clear(color='white', depth=True)
        for obj in self.objects:
            obj.sphere.draw()
        #self.sphere.draw()

if __name__ == 'main':
    win = Canvas(spindle)
    app.run()

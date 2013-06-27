# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import matplotlib.pyplot as plt


class Drawer():
    """
    """

    def __init__(self, meta_instance):
        """
        :param meta_instance: Should have been already perform a simulation
        :type meta_instance: Metaphase instance
        """

        self.meta = meta_instance
        self.KD = self.meta.KD
        self.timelapse = self.meta.timelapse
        self.num_steps = self.meta.num_steps

    def show_all(self, axes=None, fname=None):
        """
        Plot the different trajectories

        :param fname: if fname is set plot won't be shown but saved on your hard drive. Filetype is automatic according to the file extension (.png, .pdf, etc)
        :type fname: string, optional
        """
        N = int(self.KD.params['N'])
        if axes == None:
            fig = plt.figure()
            axes = fig.gca()
        spbRtraj = self.KD.spbR.traj
        spbLtraj = self.KD.spbL.traj
        axes.plot(self.timelapse, spbRtraj, color='r', ls='-', lw=1)
        axes.plot(self.timelapse, spbLtraj, color='r', ls='-', lw=1)

        fmt_list = ['g-', 'b-', 'm-']
        for n in range(N):
            ch = self.KD.chromosomes[n]
            right_traj = ch.cen_A.traj
            left_traj = ch.cen_B.traj
            fmt = fmt_list[np.mod(n, 3)]
            line1 = axes.plot(self.timelapse, right_traj, fmt, alpha=0.5)
            line2 = axes.plot(self.timelapse, left_traj, fmt, alpha=0.5)
            if n == 0:
                line1[0].set_alpha(1.)
                line2[0].set_alpha(1.)
        axes.set_xlabel('Time (seconds)', fontsize='small')
        axes.set_ylabel(u'Distance from center (um)', fontsize='small')

        if fname:
            plt.savefig(fname)
        else:
            plt.show()

    def show_one(self, n=0, fig=None, fname=None):
        """
        Shows chromosome n trajectory and plug state

        :param fname: if fname is set plot won't be shown but saved on your hard drive. Filetype is automatic according to the file extension (.png, .pdf, etc)
        :type fname: string, optional
        """
        dt = self.KD.params['dt']
        ch = self.KD.chromosomes[n]

        if fig == None:
            fig = plt.figure()
        fig.clear()

        #fig.add_subplot(312)
        gridspec = plt.GridSpec(5, 1)
        subplotspec = gridspec.new_subplotspec((1, 0), rowspan=3)
        traj_ax = fig.add_subplot(subplotspec)
        traj_ax.plot(self.timelapse, ch.cen_A.traj, 'g', lw=2, alpha=0.5)
        traj_ax.plot(self.timelapse, ch.cen_B.traj, 'purple', lw=2, alpha=0.5)
        traj_ax.plot(self.timelapse, self.KD.spbR.traj, 'k')
        traj_ax.plot(self.timelapse, self.KD.spbL.traj, 'k')
        for plugsite in ch.cen_A.plugsites:
            traj_ax.plot(self.timelapse, plugsite.traj, 'g')
        for plugsite in ch.cen_B.plugsites:
            traj_ax.plot(self.timelapse, plugsite.traj, 'purple')
        traj_ax.set_xticks([], '')

        erroneous_hist = ch.erroneous_history
        correct_hist = ch.correct_history

        subplotspec = gridspec.new_subplotspec((0, 0), rowspan=1)
        ax = fig.add_subplot(subplotspec, sharex=traj_ax)
        ax.plot(self.timelapse, erroneous_hist[:, 0], 'r',
                label='number of erroneoustellic MTs')
        ax.plot(self.timelapse, correct_hist[:, 0], 'g',
                label='number of correct MTs')
        ax.axis((0, self.num_steps * dt, -0.5, 4.5))
        ax.set_xticks([], '')

        subplotspec = gridspec.new_subplotspec((4, 0), rowspan=1)
        ax = fig.add_subplot(subplotspec, sharex=traj_ax)
        ax.plot(self.timelapse, erroneous_hist[:, 1], 'r',
                label='number of erroneous MTs')
        ax.plot(self.timelapse, correct_hist[:, 1], 'purple',
                label='number of correct MTs')
        ax.axis((0, self.num_steps * dt, -0.5, 4.5))

        if fname:
            plt.savefig(fname)
        else:
            plt.show()

# -*- coding: utf-8 -*-

"""
Usefull functions to plot with matplotlib
"""

import gc
import matplotlib.pyplot as plt


def dic_plot(plot_data, fname=None):
    """
    Plot curve with dictionnary description
    """

    d = plot_data

    if "params_box" in d.keys():
        fig = plt.figure(figsize=(12, 4))
        grid_size = (1, 6)
        ax = plt.subplot2grid(grid_size, (0, 0),
                              rowspan=1, colspan=4)
        box = plt.subplot2grid(grid_size, (0, 4),
                               rowspan=1, colspan=2)

        box.set_axis_off()
    else:
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)

    # Do we draw legend ?
    draw_legend = True
    if 'legend' in d.keys() and d['legend'] == False:
        draw_legend = False

    # Do we draw error ?
    draw_erro = True
    if 'error' in d.keys() and d['error'] == False:
        draw_erro = False

    # Add label and title
    if 'label' in d['xaxis'].keys():
        ax.set_xlabel(d['xaxis']['label'])
    if 'label' in d['yaxis'].keys():
        ax.set_ylabel(d['yaxis']['label'])
    if 'title' in d.keys():
        ax.set_title(d['title'], fontsize=16)

    # Set plot configuration
    plotter = ax.plot
    if 'logx' in d.keys() and d['logx'] == True and\
        'logy' in d.keys() and d['logy'] == True:
        plotter = ax.loglog
    if 'logx' in d.keys() and d['logx'] == True:
        plotter = ax.semilogx
    if 'logy' in d.keys() and d['logy'] == True:
        plotter = ax.semilogy

    xaxis = d['xaxis']['data']
    legends = []
    legends_label = []

    # Find axis limit
    xmin = int(d['xaxis']['data'].min())
    xmax = int(d['xaxis']['data'].max())
    ymin = 0
    ymax = 1

    for yaxis in d['yaxis']['axis']:
        mu = yaxis['data']
        color = yaxis['color']

        if 'plot_args' in yaxis.keys():
            yplot, = plotter(xaxis, mu, color=color, **yaxis['plot_args'])
        else:
            yplot, = plotter(xaxis, mu, color=color)

        if 'legend' in yaxis.keys() and draw_legend:
            legends.append(yplot)
            legends_label.append(yaxis['legend'])

        if 'error' in yaxis.keys() and draw_erro:
            sigma = yaxis['error']
            ax.fill_between(xaxis, mu + sigma, mu - sigma,
                facecolor=color, alpha=0.5)

        # Save ymin and ymax
        if mu.min() < ymin:
            ymin = mu.min()
        if mu.max() > ymax:
            ymax = mu.max()

    # Move legend outside plot
    if draw_legend:
        lgd = ax.legend(tuple(legends), tuple(legends_label),
                        loc='upper center', bbox_to_anchor=(0.5, -0.1))
    else:
        lgd = None

    # Display grid on plot
    ax.grid()

    # Set axis limit
    xmin, xmax, ymin, ymax = ax.axis()

    if 'limit_y_min' in d.keys():
        ymin = d["limit_y_min"]

    ax.axis((xmin, xmax, ymin, ymax))

    if "params_box" in d.keys():
        txtstr = "Simulation parameters\n\n"

        for p in plot_data["params_box"]:
            line = "%s : %s" % (p['name'], p['data'])
            txtstr += line + "\n"

        # these are matplotlib.patch.Patch properies
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        box.text(0, 0.9, txtstr,
                 horizontalalignment='left',
                 verticalalignment='center',
                 fontsize=12,
                 bbox=props)

    # Add annotations
    if "annotations" in d.keys():
        for ann in d["annotations"]:
            ax.annotate(**ann)

    fig.tight_layout(pad=1)

    if fname:
        if lgd:
            plt.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
        else:
            plt.savefig(fname, bbox_inches='tight')
    else:
        plt.show()

    # Cleanup memory
    fig.clf()
    plt.close()
    gc.collect()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 06-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""

import numpy as np
import matplotlib.pyplot as plt

from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mpy

import h5py


def make_movie(h5filename, fps=20, output='movie.mp4'):
    with h5py.File(h5filename) as f:
        t = f['1/t'][:]
        x = f['1/x'][:]
        psiI = f['1/psiI'][:]
        psiR = f['1/psiR'][:]

    n = np.hypot(psiR, psiI)
    margin = 0.05
    ylim = n.min() - margin * n.ptp(), n.max() + margin * n.ptp()
    print(t.shape, x.shape, n.shape)
    Nframes = len(t)
    duration = Nframes / fps

    fig, ax = plt.subplots()
    line, = ax.plot(x, n[0])
    ax.plot(x, n[-1], ls='--', alpha=0.6)
    ax.set_ylim(ylim)
    plt.show()

    def make_frame_mpl(_t):
        ix = int(_t/duration * Nframes)
        # print(ix)
        line.set_ydata(n[ix])  # <= Update the curve
        return mplfig_to_npimage(fig)  # RGB image of the figure

    animation = mpy.VideoClip(make_frame_mpl, duration=duration)
    animation.write_videofile(output, fps=fps)
    # animation.write_gif("movie.gif", fps=20)

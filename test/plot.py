import matplotlib.pyplot as plotpy
import numpy


def spectrogram(dfts, sr, hopsize=1, xlim=None, ylim=None, clim=-120, cmap='inferno', ylog=False):

    def lim():

        if xlim is not None:
            if isinstance(xlim, (list, tuple)):
                plotpy.xlim(xlim)
            else:
                plotpy.xlim(0, xlim)

        if ylim is not None:
            if isinstance(ylim, (list, tuple)):
                plotpy.ylim(ylim)
            else:
                plotpy.ylim(0, ylim)

        if clim is not None:
            if isinstance(clim, (list, tuple)):
                plotpy.clim(clim)
            else:
                plotpy.clim(clim, 0)

    with numpy.errstate(divide='ignore', invalid='ignore'):
        dfts = 20 * numpy.log10(numpy.abs(dfts))

    time = numpy.array([hop * hopsize / sr for hop in range(dfts.shape[0])])
    freq = numpy.linspace(0, sr / 2, dfts.shape[1])

    extent = (numpy.min(time), numpy.max(time), numpy.min(freq), numpy.max(freq))

    plotpy.imshow(dfts.T, aspect='auto', cmap=cmap, extent=extent, interpolation='nearest', origin='lower')
    colorbar = plotpy.colorbar()

    plotpy.xlabel('s')
    plotpy.ylabel('Hz')
    colorbar.set_label('dB')

    if ylog:
        plotpy.yscale('log')

    lim()

    return plotpy

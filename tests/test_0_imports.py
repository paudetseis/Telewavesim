def test_numpy_import():
    import numpy


def test_obspy_import():
    import obspy


def test_telewavesim_modules():
    import telewavesim
    from telewavesim import elast
    from telewavesim import utils
    import matplotlib
    matplotlib.use('Agg')
    from telewavesim import wiggle
    from telewavesim.rmat_f import conf
    from telewavesim.rmat_f import plane

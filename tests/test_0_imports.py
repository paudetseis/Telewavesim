import telewavesim

## ===================


def test_numpy_import():
    import numpy

    return


def test_obspy_import():
    import obspy

    return


def test_pyfftw_import():
    import pyfftw

    return


def test_telewavesim_modules():
    import telewavesim
    from telewavesim import conf
    from telewavesim import elast
    from telewavesim import utils
    from telewavesim import wiggle
    from telewavesim.rmat_f import conf
    from telewavesim.rmat_f import plane


if __name__ == "__main__":
    test_numpy_import()
    test_obspy_import()
    test_pyfftw_import()
    test_telewavesim_modules()
#!/usr/bin/env python

import glob

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('plenslike',parent_package,top_path)
    config.add_extension('_plenslike', ['plenslike/plenslike.c', 'plenslike/wignerd.c', 'plenslike/plenslike.pyf'], undef_macros=['NDEBUG'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(name='plenslike',
          packages=['plenslike'],
          package_data={'plenslike':
                        [tf.replace("plenslike/", "") for tf in glob.glob('plenslike/data/v50/plenslike_v50_mono_???*.dat')] },
          configuration=configuration)
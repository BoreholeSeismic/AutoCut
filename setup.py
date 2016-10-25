#!/usr/bin/env python

"""setup.py: setuptools control."""

import re
from setuptools import setup

version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('autocut/auto_cut.py').read(),
    re.M
    ).group(1)

setup(name='cmd-autocut',
      packages=['autocut'],
      entry_points = {
        "console_scripts": ['autocut = autocut.auto_cut:main']
        },
      version=version,
      description='Slice long seismic traces into short events',
      install_requires=[],
      url='https://github.com/siruix/autocut',
      author='Sirui Xing',
      author_email='sirui.xing@gmail.com',

     )
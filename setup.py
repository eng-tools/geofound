from setuptools import setup, find_packages

setup(name='geofound',
      version='0.2.0',
      description='A package to assess the bearing capacity and settlement of geofound',
      url='',
      author='Maxim Millen',
      author_email='mmi46@uclive.ac.nz',
      keywords=['geotechnical', 'engineering', 'structural'],
      license='MIT',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Scientific community',
          'Topic :: Civil Engineering',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
      ],
      packages=find_packages(exclude=['contrib', 'docs', 'tests']),
      install_requires=[],
      # List additional groups of dependencies here (e.g. development
      # dependencies). You can install these using the following syntax,
      # for example:
      # $ pip install -e .[dev,test]
      extras_require={
          'test': ['pytest'],
      },
      python_requires='>=3',
      package_data={
          'models': ['models_data.dat'],
      },
      zip_safe=False)


# From python packaging guides
# versioning is a 3-part MAJOR.MINOR.MAINTENANCE numbering scheme,
# where the project author increments:

# MAJOR version when they make incompatible API changes,
# MINOR version when they add functionality in a backwards-compatible manner, and
# MAINTENANCE version when they make backwards-compatible bug fixes.

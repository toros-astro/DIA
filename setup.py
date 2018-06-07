from setuptools import setup, Extension
import numpy
import os

# Get the version from astroalign file itself (not imported)
with open('oisdiff.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            _, _, module_version = line.replace("'", '').split()
            break

oisdiffmodule = Extension('oismodule', 
                    sources=[os.path.join('src', f) for f in ('oismodule.c', 'oisdifference.c')],
                    extra_compile_args=["-std=c99"],
                    libraries = ['m',],
)

setup(name='oisdiff',
      version=module_version,
      description='Optimal Image Subtraction',
      author='Martin Beroiz & Ryan Oelkers',
      author_email='martinberoiz@gmail.com',
      url='https://github.com/torosastro/DIA',
      py_modules=['oisdiff', ],
      ext_modules=[oisdiffmodule],
      include_dirs=[numpy.get_include()],
      install_requires=["numpy>=1.6"],
      test_suite='tests',
      )

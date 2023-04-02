from setuptools import Extension, setup, find_packages

module = Extension("mykmeanssp", sources=['spkmeansmodule.c', 'spkmeans.c','main.c'])
setup(name='mykmeanssp',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])
#!/usr/bin/env python
import sys
from os import path, name as osname
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

sources = ['simplesum.c', \
           'rotatesum.c', \
           'fastincommsum.c', \
           'vec3.c', \
           'mat3.c', \
           'pile.c', \
           'dipolartensor.c',\
           'utility.c']
src_sources = []
for s in sources:
    src_sources.append(
    path.join('.','src',s)
    )


openmp_compile_args = []
openmp_link_args = []


## In case of missing numpy array headers
try:
    from numpy import get_include as numpy_get_include
    numpy_include_dir = [numpy_get_include()]
except:
    numpy_include_dir = []

# Ugly hack to set compiler flags 
COMPILE_ARGS = {'msvc':[],'gcc':[],'unix':[]}
LINK_ARGS = {'msvc':[],'gcc':[],'unix':[]}

for compiler, args in [
        ('msvc', ['/EHsc', '/DHUNSPELL_STATIC']),
        ('gcc', ['-O3', '-g0', '-std=c89']),
        ('unix', ['-O3', '-g0', '-std=c89'])]:
    COMPILE_ARGS[compiler] += args

# add m lib only on linux/osx
for compiler, args in [
        ('msvc', []),
        ('unix', ['-lm']),
        ('gcc', ['-lm'])]:
    LINK_ARGS[compiler] += args  

# Ugly hack to have openMP as option
if "--with-openmp" in sys.argv:
    for compiler, args in [
            ('msvc', ['/openmp']),
            ('unix', ['-fopenmp']),
            ('gcc', ['-fopenmp'])]:
        COMPILE_ARGS[compiler] += args    
    for compiler, args in [
            ('msvc', []),
            ('unix', ['-lgomp']),
            ('gcc', ['-lgomp'])]:
        LINK_ARGS[compiler] += args    

    sys.argv.remove("--with-openmp")

# Ugly hack to have direct library link as option
for carg in sys.argv:
    if "--with-lfc-lib" in carg:
        if osname != 'posix':
            print("Static/dynamic lib linking not supported.")
            break
        
        src_sources = []
            
        for compiler, args in [
                ('msvc', []),
                ('unix', [carg.split('=')[1]]),
                ('gcc', [carg.split('=')[1]])]:
            LINK_ARGS[compiler] += args    
            LINK_ARGS[compiler] += "-lgsl -lgslcblas".split()
    
        sys.argv.remove(carg)
        break

class build_ext_compiler_check(build_ext):
    def build_extensions(self):
        
        compiler = self.compiler.compiler_type
        cargs = COMPILE_ARGS[compiler]
        for ext in self.extensions:
            ext.extra_compile_args = cargs
            
        largs = LINK_ARGS[compiler]
        for ext in self.extensions:
            ext.extra_link_args = largs
        
        build_ext.build_extensions(self)



    
setup(name='LFC',
      version='0.1',
      description='Local Field Components (or lighting fast calculator) for muesr package',
      author='Pietro Bonfa',
      author_email='pietro.bonfa@fis.unipr.it',
      url='https://github.com/bonfus/muesr',
      packages=['LFC',],
      ext_modules=[Extension('lfclib', sources = ['python/lfclib.c',]+src_sources,
                                      libraries=[],
                                      include_dirs=numpy_include_dir)],
     package_dir={'LFC': 'python' },
     install_requires=[
          'numpy >= 1.6',
     ],
     test_suite="tests",
     cmdclass={ 'build_ext': build_ext_compiler_check }
     )

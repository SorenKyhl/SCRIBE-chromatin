from setuptools import find_packages, setup, Extension
import subprocess
import platform
from pybind11.setup_helpers import Pybind11Extension

# Extra compile args for C++17
extra_compile_args = ['-std=c++17']
extra_link_args = []

# macOS-specific linker flags
if platform.system() == 'Darwin':
    extra_link_args.append('-undefined')
    extra_link_args.append('dynamic_lookup')

module1 = Pybind11Extension(
    name='pylib.pyticg',
    include_dirs=['include'],
    language='c++',
    sources=["src/pybind_Sim.cpp"],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
)

setup(
    name='pylib',
    packages=find_packages(include=['pylib', 'pylib.*']),
    package_data={
        'pylib': ['defaults/*.json'],
    },
    include_package_data=True,
    version='0.1.3',
    description='SCRIBE Python library',
    author='Soren Kyhl',
    license='MIT',
    install_requires=['hic-straw', 'jsbeautifier'],
    ext_modules=[module1],
)

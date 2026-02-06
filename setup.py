import platform
import subprocess
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension
from setuptools import find_packages, setup


# Write git commit hash into scribe/_version_info.py so it's available at runtime
def _write_version_info():
    try:
        commit = (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
        dirty = (
            subprocess.call(
                ["git", "diff", "--quiet"],
                stderr=subprocess.DEVNULL,
            )
            != 0
        )
        if dirty:
            commit += "-dirty"
    except Exception:
        commit = "unknown"
    path = Path(__file__).parent / "scribe" / "_version_info.py"
    path.write_text(f'__git_commit__ = "{commit}"\n')


_write_version_info()

# Extra compile args for C++17
extra_compile_args = ["-std=c++17", "-Wno-deprecated-declarations"]
extra_link_args = []

# macOS-specific linker flags
if platform.system() == "Darwin":
    extra_link_args.append("-undefined")
    extra_link_args.append("dynamic_lookup")

module1 = Pybind11Extension(
    name="scribe.scribe_engine",
    include_dirs=["include"],
    language="c++",
    sources=["src/pybind_Sim.cpp"],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
)

setup(
    name="scribe",
    packages=find_packages(include=["scribe", "scribe.*"]),
    package_data={
        "scribe": ["defaults/*.json"],
    },
    include_package_data=True,
    version="0.1.3",
    description="SCRIBE Python library",
    author="Soren Kyhl",
    license="MIT",
    install_requires=["hic-straw", "jsbeautifier"],
    ext_modules=[module1],
)

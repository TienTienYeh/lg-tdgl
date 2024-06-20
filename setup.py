"""
# LG-TDGL


"""

from setuptools import find_packages, setup

DESCRIPTION = "LGTDGL: ????? Time-dependent Ginzburg-Landau in Python."
LONG_DESCRIPTION = __doc__

NAME = "lgtdgl"
AUTHOR = "TTY"
AUTHOR_EMAIL = "ttyeh0212@gmail.com"
URL = "https://github.com/TienTienYeh/lg-tdgl"
LICENSE = "MIT"
PYTHON_VERSION = ">=3.8, <3.12"

INSTALL_REQUIRES = [
    "cloudpickle",
    "h5py",
    "joblib",
    "jupyter",
    "matplotlib",
    "meshpy",
    "numba",
    "numpy",
    "pint==0.23",
    "pytest",
    "pytest-cov",
    "scipy<1.11",
    "shapely",
    "tqdm",
    "tdgl==0.8.2",
]

EXTRAS_REQUIRE = {
    "dev": [
        "black",
        "isort",
        "pre-commit",
    ],
    "docs": [
        "IPython",
        # https://github.com/readthedocs/sphinx_rtd_theme/issues/1115
        "sphinx==5.3.0",
        "sphinx-rtd-theme>=0.5.2",
        "sphinx-autodoc-typehints",
        "nbsphinx",
        "pillow",
        "sphinx_toolbox",
        "enum_tools",
        "sphinx-argparse",
        "sphinxcontrib-bibtex",
    ],
    "umfpack": [
        "swig",
        "scikit-umfpack",
    ],
    "pardiso": [
        "pypardiso",
    ],
}

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: MacOS",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: Microsoft :: Windows",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Physics",
]

PLATFORMS = ["Linux", "Mac OSX", "Unix", "Windows"]
KEYWORDS = "superconductor vortex Ginzburg-Landau Laguerre Gaussian Light"


setup(
    name=NAME,
    version=__version__,  # noqa: F821
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    license=LICENSE,
    packages=find_packages(),
    include_package_data=True,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    keywords=KEYWORDS,
    classifiers=CLASSIFIERS,
    platforms=PLATFORMS,
    python_requires=PYTHON_VERSION,
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
)

from setuptools import setup, find_packages, Extension

import numpy as np

extensions = [
    Extension(
        "spiir.search.skymap._skymap",
        sources=["src/spiir/search/skymap/_skymap.c"],
        include_dirs=[np.get_include()],
        libraries=["lal", "gsl"],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
    ),
]


setup(
    name="spiir",
    version="0.0.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "wheel",
        "setuptools",
        "lalsuite",
        "ligo.skymap",
        "astropy",
        "astropy-healpix",
        "python-ligo-lw==1.8.1",
        "igwn-alert",
        "ligo-gracedb",
        "toml",
        "numpy",
        "scipy",
        "pandas",
        "matplotlib",
    ],
    extras_require={
        "tensorflow": ["tensorflow>=2.8", "tensorflow-probability", "scikit-learn"],
        "torch": ["torch", "torchaudio", "scikit-learn"],
        "pycbc": ["pycbc"],
    },
    ext_modules=extensions,
    include_package_data=True,
    description="A Python library for the SPIIR gravitational wave science pipeline.",
    author="Daniel Tang",
    author_email="daniel.tang@uwa.edu.au",
)

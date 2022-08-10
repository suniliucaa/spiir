from setuptools import setup, find_packages

setup(
    name="spiir",
    version="0.0.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "wheel",
        "setuptools",
        "lalsuite",
        "astropy",
        "gwpy",
        "python-ligo-lw",
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
    description="A Python library for the SPIIR gravitational wave science pipeline.",
    author="Daniel Tang",
    author_email="daniel.tang@uwa.edu.au",
)

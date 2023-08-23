from setuptools import setup, find_packages

VERSION = "0.0.1"
DESCRIPTION = "BeyondGR ID solver."
LONG_DESCRIPTION = (
    "Collection of python scripts for initial data solves in beyond GR theories."
)

setup(
    name="BeyondGRID",
    version=VERSION,
    author="Peter James Nee",
    author_email="peter.nee@aei.mpg.de",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=["numpy"],
    keywords=["python", "GR", "sGB"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
    ],
)

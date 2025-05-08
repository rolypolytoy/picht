from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="picht",
    version="1.1.4",
    author="Rishiit Sharma",
    author_email="rishiitsharma@gmail.com",
    description="Electron and ion optics simulation using the Finite Difference Method (FDM)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rolypolytoy/picht",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Intended Audience :: Science/Research",
        "Development Status :: 5 - Production/Stable",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy==1.26.4",
        "matplotlib",
        "scipy",
        "numba",
        "mendeleev",
    ],
    license="MIT",
)
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyProfiler",
    version="0.0.1",
    author="Laurent K, David M",
    author_email="dmoi@unil.ch",
    description="Phylogenetic profiling package",
    package_dir={'':'.'},
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: MIT License",
        "Operating System :: OS Independent",
    ],
)

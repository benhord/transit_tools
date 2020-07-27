import setuptools

with open("README.md", "r") as fh:
    long_desccription = fh.read()

setuptools.setup(
    name="transit_tools_benhord",
    version="0.0.1",
    author="Ben Hord",
    author_email="benhord@astro.umd.edu",
    description="Some exoplanet transit tools I frequently use.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/benhord/transit_tools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT Licencse",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)

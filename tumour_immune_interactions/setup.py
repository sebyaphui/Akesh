from setuptools import setup, find_packages  # noqa: D100

setup(
    name="tumour_immune_interactions",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
    ],
)

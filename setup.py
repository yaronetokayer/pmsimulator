import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pmsimulator',
    version='0.0.1',
    author='Yarone Tokayer',
    author_email='yarone.tokayer@yale.edu',
    description='A gravitational N body code for Python using a particle mesh method',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/yaronetokayer/pmsimulator',
    project_urls = {
        "Bug Tracker": "https://github.com/yaronetokayer/pmsimulator/issues"
    },
    license='GNU GPLv3',
    packages=[
        'pmsimulator'
    ],
    install_requires=[
        'numpy',
        'matplotlib'
    ],
)

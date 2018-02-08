from setuptools import setup, find_packages
import kamikaze

# Add setuptools boilerplate
setup(
    name='kamikaze',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=kamikaze.__version__,

    description='Oligo design tool for kf system',
    long_description='design tool for kf system',

    # The project's main homepage.
    url='https://github.com/elijahc/kamikaze',

    # Author details
    author='Elijah Christensen',
    author_email='ejd.christensen@gmail.com',

    # Choose your license
    license='Apache 2.0',
    packages=find_packages(),
    package_data={'mypkg': ['data/*.fa']},
)


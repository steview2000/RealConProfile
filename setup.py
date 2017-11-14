try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md', 'r') as f:
    long_description = f.read()

setup(name='RealConProfileSingle',
      py_modules=['RealConProfileSingle'],
      version='0.1',
      license='MIT',
      description='Calculates temperature profile and conductive heat transport for RBC',
      long_description=long_description,
      author='Stephan Weiss ',
      author_email='stephan.weiss@ds.mpg.de',
      url='',
      download_url='',
      keywords=['NIST','REFPROP','Fluid properties','convection'],
      install_requires=['numpy'],
      classifiers=["Topic :: Scientific/Engineering",
                   "License :: OSI Approved :: MIT License",
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3.5",
                   "Development Status :: 4 - Beta"],
      )

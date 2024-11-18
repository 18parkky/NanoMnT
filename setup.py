from setuptools import setup

setup(
    name='NanoMnT',
    version='0.1.0',    
    description='A collection of Python scripts for analyzing STR from ONT data',
    url='https://github.com/18parkky/NanoMnT',
    author='Gyumin Park',
    author_email='18parkky@gm.gist.ac.kr',
    license='BSD 2-clause',
    packages=['NanoMnT'],
    entry_points={
        'console_scripts' : [ 'getAlleleTable = NanoMnT.getAlleleTable:main', 
                             'getLocusTable = NanoMnT.getLocusTable:main', 
                             'findInformativeLoci = NanoMnT.findInformativeLoci:main', 
                             ] 
    },
    install_requires=['pysam>=0.20.0',
                      'Levenshtein>=0.21.1',     
                      'numpy==1.26.4',
                      'pandas>=2.0.0',
                      'scipy>=1.7.1',
                      'matplotlib>=3.7.1',
                      'seaborn>=0.13.0',           
                      ],

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3.9.7',
    ],
)
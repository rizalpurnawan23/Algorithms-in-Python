from setuptools import setup, find_packages

# Package metadata
NAME = 'Algorithms-in-Python'
VERSION = '1.0'
AUTHOR = 'Rizal Purnawan'
DESCRIPTION = 'A Python package for multiple purposes'
URL = 'https://github.com/rizalpurnawan23/Algorithms-in-Python'
LICENSE = 'MIT'
CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Practitioners',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
]
KEYWORDS = 'python package example'
INSTALL_REQUIRES = [
    'numpy>=1.0',        # Requires NumPy version 1.0 or higher
    'pandas>=1.0',       # Requires Pandas version 1.0 or higher
    'matplotlib>=3.0',   # Requires Matplotlib version 3.0 or higher
    'seaborn>=0.10',     # Requires Seaborn version 0.10 or higher
    'scikit-learn>=0.24',# Requires Scikit-learn version 0.24 or higher
    'tensorflow>=2.0',   # Requires TensorFlow version 2.0 or higher
]


# Load long description from README file
# with open('README.md', 'r', encoding='utf-8') as f:
#     LONG_DESCRIPTION = f.read()

# Setup configuration
setup(
    name= NAME,
    version= VERSION,
    author= AUTHOR,
    description= DESCRIPTION,
    # long_description=LONG_DESCRIPTION,
    # long_description_content_type= 'text/markdown',
    url= URL,
    license= LICENSE,
    classifiers= CLASSIFIERS,
    keywords= KEYWORDS,
    packages= find_packages(),
    install_requires= INSTALL_REQUIRES,
    python_requires= '>=3.7, <4',
)

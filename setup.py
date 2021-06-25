from setuptools import setup, find_packages

setup(
    name='scanit',
    version='0.1',
    packages=find_packages(exclude=['tests*']),
    description='Representation learning for spatial transcriptomics data.',
    author='Zixuan Cang',
    author_email='cangzx@gmail.com'
)

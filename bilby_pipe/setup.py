from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name="testpackage",
     version="1.0",
     description="test for sine gaussian waveform",
     author="Hemanta_ph",
     packages=['testpackage'],
     setup_requires=["numpy"],
     install_requires=['numpy'],
     )
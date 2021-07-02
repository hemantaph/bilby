from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()
    
setup(name='mydemo',
     version='0.0.1',
     description='Demo',
     long_description=readme(),
     long_description_content_type='text/markdown',
     classifiers=[
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Operating System :: OS Independent"],
      url='',
      author='Hemanta',
      author_email='hemantaphurailatpam@gmail.com',
      keywords='core package',
      licence='CUHK',
      packages=['mydemo','mydemo.addition'],
      install_requires=['cython','numpy','bilby','gwsurrogate'],
      include_package_data=True,
      zip_safe=False)
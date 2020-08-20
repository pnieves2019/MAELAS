import re
from setuptools import setup
from os import makedirs,environ
from subprocess import call
import shlex
from sys import argv,path,version_info

packages_maelas = ['maelas']
executables_maelas = ['maelas/bin/maelas']

requirements_maelas = []

VERSION='1.0.2'

if __name__ == '__main__':
    myself=environ['PWD']
    options=""
    try:
        if argv[1] != 'install':
            exit(-1)
    except IndexError:
        exit(-1)
    for arg in argv[2:]:
        options+=" %s"%arg

    options=re.sub('=',' ',options)
    options=re.sub('~',environ['HOME'],options)
    if '--install_reqs' in argv:
        print('Installing requirements')
        call(shlex.split('./install-requirements.sh'+options),shell=False)
        argv.remove('--install_reqs')

    if '--user' in argv:
        try:
            locDir = environ['PYTHONUSERBASE']
            with open(environ['HOME']+'/.bashrc','a+') as bashrc:
                bashrc.write('\nexport PYTHONPATH=%s/lib/python%d.%d/site-packages/maelas-%s-py%d.%d.egg:$PYTHONPATH\n'%(locDir,version_info[0],version_info[1],VERSION,version_info[0],version_info[1]))
            makedirs(locDir,exist_ok=True)
        except KeyError:
            locDir = environ['HOME']+'/.local'
            print('Installing in %s'%(environ['HOME']))
        options+=' --prefix %s'%locDir
    try:
        with open('requirements.txt','r') as reqs:
            for line in reqs.readlines():
                package = re.sub('[><=]+.*','',line)
                package = re.sub('\s','',package)
                if len(package)>0:
                    requirements_maelas.append(package)
    except FileNotFoundError:
        requirements_maelas = [ "docutils", "pymatgen", "scikit-learn", "pyfiglet", "argparse", "numpy", "matplotlib", "scipy", "setuptools" ]

    setup(name='maelas',
          version=VERSION,
          description='This is a software to calculate anisotropic magnetostriction coefficients and magnetoelastic constants',
          long_description="""
          This is a software to calculate anisotropic magnetostriction coefficients and magnetoelastic constants
          """,
          author='P. Nieves, S. Arapan, S.H. Zhang, A.P. Kądzielawa, R.F. Zhang and D. Legut',
          author_email='pablo.nieves.cordones@vsb.cz',
          keywords="Magnetostriction, Magnetoelasticity, High-throughput computation, First-principles calculation",
          url="https://github.com/pnieves2019/MAELAS",   # project home page, if any
          project_urls={
              "Documentation": "https://github.com/pnieves2019/MAELAS/Manual.pdf",
              "Source Code": "https://github.com/pnieves2019/MAELAS/maelas.py",
              "Examples": "https://github.com/pnieves2019/MAELAS/Examples",
          },
          classifiers=[
              "License :: OSI Approved :: BSD License"
          ],
          python_requires='>=3.6',   
          packages=packages_maelas,
          install_requires=requirements_maelas,
          provides=['maelas'],
          scripts=executables_maelas,
          platforms=['POSIX'])







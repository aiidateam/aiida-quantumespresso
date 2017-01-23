from setuptools import setup, find_packages

if __name__ == '__main__':
    setup(
        version='0.7.1b0',
        name='aiida_quantumespresso',
        url='http://www.aiida.net',
        license='MIT License',
        author='The AiiDA team',
        author_email='developers@aiida.net',
        include_package_data=True,
        classifiers=[
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 2.7',
            'Development Status :: 4 - Beta',
        ],
        install_requires=[
            'aiida[ssh]'
        ],
        packages=find_packages(),
        entry_points={
            'aiida.calculations': [
                'quantumespresso.pw = aiida_quantumespresso.calculations.pw:PwCalculation',
                'quantumespresso.cp = aiida_quantumespresso.calculations.cp:CpCalculation',
                'quantumespresso.pwimmigrant = aiida_quantumespresso.calculations:PwimmigrantCalculation',
            ],
            'aiida.parsers': [
              'quantumespresso.basicpw = aiida_quantumespresso.basicpw:BasicpwParser',
              'quantumespresso.cp = aiida_quantumespresso.plugins.cp:CpParser',
              'quantumespresso.pw = aiida_quantumespresso.plugins.pw:PwParser',
            ],
            'aiida.tools.dbexporters.tcod_plugins': [
              'cp = aiida_quantumespresso.dbexport.tcodplugins.cp:CpTcodtranslator',
              'pw = aiida_quantumespresso.dbexport.tcodplugins.pw:PwTcodtranslator'
            ]
        }
    )

==========
mpcorbfile
==========


.. image:: https://img.shields.io/pypi/v/mpcorbfile.svg
        :target: https://pypi.python.org/pypi/mpcorbfile


.. image:: https://readthedocs.org/projects/mpcorbfile/badge/?version=latest
        :target: https://mpcorbfile.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

Read and write MPCORB files from  minor planet center.

Features
--------
This package provides a simple interface to read and write MPCORB files from the minor planet center formated 
as described in https://minorplanetcenter.net/iau/info/MPOrbitFormat.html

Its not provides any kind of orbital calculation, just read and write the files. See Usage for more details.

Reading::

        import mpcorbfile
        
        #read MPCORB.DAT download  from:
        #https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT
        mpc = mpcorbfile.mpcorb_file('MPCORB.DAT')
        #mpc.bodies it's a list of dictionaries with the data
        print(f'{mpc.bodies[0]})
        {'packed_designation': 'a7943',
        'H': 24.19,
        'G': 0.15,
        'Epoch': datetime.datetime(2025, 5, 5, 0, 0),
        'M': 184.7538,
        'Peri': 195.63684,
        'Node': 146.91033,
        'i': 11.6102,
        'e': 0.0894684,                 
        'n': 1.13494811,
        'a': 0.9102319,
        'U': '2',
        'Ref': 'E2024-JU2',
        'Num_obs': 1055,
        'Num_opps': 2,
        'Arc_length': '2012-2013',
        'rms': 0.36,
        'Perturbers': 'M-v',
        'Perturbers_2': '3Ek',
        'Computer': 'Veres',
        'Hex_flags': 2049,
        'Number': '(367943)',
        'Name': 'Duende',
        'Last_obs': '20130221',
        'epochJD': 2460800.5,
        'designation': '367943',
        'discover_date': nan,
        'orbit_type': 'NEO;Athen'}

* Free software: MIT license
* Documentation: https://mpcorbfile.readthedocs.io.


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

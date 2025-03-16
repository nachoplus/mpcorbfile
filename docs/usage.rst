=====
Usage
=====

CLI commands:

* mpcorbfile2json: reads an MPCORB file and writes a json file with the same data
* json2mpcorb: reads a json file and writes an MPCORB file with the same data


API interface to read and write MPCORB files from python code examples

Reading::

        import mpcorbfile
        
        #read MPCORB.DAT download  from:
        #https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT
        mpc = mpcorbfile.mpcorb_file('MPCORB.DAT')

        
        #alternatively you can use the mpcorb_extended.json.gz
        #download from: 
        #https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz
        mpc = mpcorbfile.mpcorb_file()
        bodies=mpc.read_json('mpcorb_extended.json') # mpc.bodies contains the data but also returns it

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

Writing::

        import mpcorbfile

        import datetime
        mpc=mpcorbfile.mpcorb_file()
        Name='2026 AA1'
        epoch = datetime.datetime(2026,1,1)
        packed_designation=mpc._pack_designation(Name)
        body={'packed_designation':packed_designation,'H': 3.2, 'G': 0.15,
                'Epoch': epoch, 'M': 0.0, 'Peri': 0.0, 'Node': 0.0, 
                'i': 0.0, 'e': 0.0, 'n': 0.0, 'a': 0.0,'U': 0.0,
                'Ref': 'MPC 12345', 'Num': 12345, 'Name': 'Test'}
        mpc.add(body)
        #write MPCORB.DAT format
        mpc.write('test.dat')
        #or write json
        mpc.write_json('test.json')



See jupyter notebooks for more examples

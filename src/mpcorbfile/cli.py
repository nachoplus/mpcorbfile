"""Console script for mpcorbfile."""
import mpcorbfile
import click as clk

    
@clk.command()
@clk.argument('mpcorbfile_name')
@clk.argument('output')
def mpcorbfilecli(mpcorbfile_name:str,output:str)->bool:
    """Console script for mpcorbfile."""
    f=mpcorbfile.mpcorb_file()
    f.read(mpcorbfile_name)
    f.write_json(output)
    return True


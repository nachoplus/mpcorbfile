"""Console script for mpcorbfile."""

import mpcorbfile
import click as clk


@clk.command()
@clk.argument("mpcfile")
@clk.argument("jsonfile")
def mpcorbfilecli(mpcfile: str, jsonfile: str) -> bool:
    """Console script for mpcorbfile."""
    f = mpcorbfile.mpcorb_file()
    f.read(mpcfile)
    f.write_json(jsonfile)
    return True

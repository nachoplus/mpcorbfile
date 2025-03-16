"""Console script for mpcorbfile."""

import mpcorbfile
import click as clk


@clk.command()
@clk.argument("mpcfile")
@clk.argument("jsonfile")
def cli_mpcorb2json(mpcfile: str, jsonfile: str) -> bool:
    """Convert MPCORB.DAT file to JSON file"""
    f = mpcorbfile.mpcorb_file()
    f.read(mpcfile)
    f.write_json(jsonfile)
    return True


@clk.command()
@clk.argument("jsonfile")
@clk.argument("mpcfile")
def cli_json2mpcorb(jsonfile: str, mpcfile: str) -> bool:
    """Convert JSON file to MPCORB.DAT file"""
    f = mpcorbfile.mpcorb_file()
    f.read_json(jsonfile)
    f.write(mpcfile, header=f"     Converted from {jsonfile} file")
    return True

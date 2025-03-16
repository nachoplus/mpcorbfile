=======
History
=======

1.0.4 (2025-03-17)
------------------
* better Documentation
* include jupyter notebook examples


1.0.3 (2025-03-16)
------------------
* better Documentation
* first sphinx documentation
* more optional fields
* change some function from private to public and viceversa
* fix get_chunk fn

1.0.2 (2025-03-10)
------------------
* json read support
* packed_designation form designation
* cli command to convert from json to mpcorb ASCII


1.0.1 (2025-03-09)
------------------
* Read discovery date from 'name' in instead of from 'packed_designation'. Packed designation losses his date meaning when asteroid get numbered 'name' field has discovery date meaning while it is provisional (no given name).
* tqdm to show progress
* fix json date serialization
* format float,integer and hex acording specification


1.0.0 (2025-03-07)
------------------
* Change internal. All objects now are store with its type (not strings) thus to_numeric() is not used ever
* Put convinience fn outside of main class

0.3.2 (2025-03-05)
------------------

* Add calculated fields (expanded designation,orbit type, date of discovery from name)
* fn to add bodies to REBOUND simulation
* Added examples notebook for use with pyephem and REBOUND
* get_chunk fn to divide bodies in packet for multiprocessing
* Bug Fix


0.1.0 (2025-03-04)
------------------

* First release on PyPI.

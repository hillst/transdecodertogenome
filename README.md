transdecodertogenome
====================

Converts transdecoder output (Gene coordinates) to genome coordinates.

##Purpose

The tool transdecoder takes generated transcripts and identifies the coding reigons within them, however, the coordinates
of the output GFF file are no longer in genome coordinates. The tool provided by trandecoder to convert the coordinates doesn't
work, so this tool will convert them. It also performs various tests to make sure the data is valid, for instance, it verifies 
the strandedness of each transcript as marked in the GFF3.

There is a test suite associated with this application as well, with the same name, but test appended to the end.

#### Usage

`python --cigar_path <path to cigar input file> --query_path <path to query file>`

#### Sample Input/Output Files

Input:  `test_query_in.txt`

Output `query_results.txt`

#### Description

This program accepts a an input file containing a list of transcripts with their chromosome locations, starting positions in terms of the **reference sequence**, and a CIGAR string describing the characteristics of the alignment.  No sequence information is given.

The query file contains a list of transcripts and the non-zero coordinate on the **transcript**

When the 4-column query file is loaded, the information is converted into a Pandas dataframe and the mapping coordinates for each position on each transcription are determined by iterating through each CIGAR operation while adjusting for the non-zero starting reference coordinate.  Stretches of the transcript that are insertions are acconted for by adding 1 to each transcript coordinate included in the region, while regions in which there is a deletion add 1 to each reference coordinate.

A dictionary of transcript/reference coordinates are returned as a dictionary, which is queryable given a transcript name and zero-based coordinate on the transcript.  All queries are fetched automatically by running the module and are written to `query_results.txt`

This program only accepts the CIGAR operations **M, X, I,** and **D**.  Other operations are not accepted and will raise an `UnsupportedOperation ` exception.  CIGAR strings that do no adhere to the accepted format will raise a `ValueError`.

If a query transcript is not found within the mapping dictionary, it will be noted as such in the output file (see example output)

The key strength of this solution is that all mapped coordinates  for an entire transcript are calculated  and returned before querying.  Using key:value pairs, it is simple to retrieve a corresponding refrence coordinate.  In addition, if a query position is out of range (i.e. a position past the end of the aligned transcript), a `KeyError` will catch the error and report it in the final output.

The instantiation of the `QueryHandler` and `CigarMapper` classes keep the query and mapping information separate from one another.  This makes debugging easier.

Although this is a straightforward and reliable solution for a small number of transcripts, it may not be the most efficient if a very large number of CIGAR sequences need to be mapped.  Each transcript coordinate has a correponding reference coordinate stored in memory. Large dictionaries and dataframes consume memory which is not the most optimal solution when large amounts of data need to be processed.  Also, finding maps though iteration over each CIGAR string is not the fastest method.  

**Testing**

This program was tested to properly handle exceptions for unsupported and invalid CIGAR strings.  The mapping algorithm was tested by manually generating mapping dictionaries for a selection of test transcripts.  These were compared against the mapping dictionaries outputted by `CigarMapper._map_positions`


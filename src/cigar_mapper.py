import argparse
import pandas as pd
from HTSeq import parse_cigar

supported_cigar_ops = ('M', 'D', 'I', 'X')
df_colnames = ['transcript_id', 'chromosome', 'ref_pos', 'cigar_string']
out_file_name = 'query_results.txt'

class UnsupportedOperation(Exception):

    def __init__(self, message=f'Unsupported cigar operation detected. \ '
                               f'Allowed operations: {supported_cigar_ops}'):
        self.message = message
        super().__init__(self.message)


class CigarMapper:

    def __init__(self, map: pd.DataFrame):
        self.input = map
        self.position_mappings = self.get_all_positions(map)

    def check_for_unsupported_ops(self, cigar_obj: list):
        """Checks for operation characters not supported by this program."""
        ops = set([op.type for op in cigar_obj])
        if not all(op in supported_cigar_ops for op in ops):
            raise UnsupportedOperation()

    def _map_positions(self, cigar_string: str, ref_position: int) -> dict:
        """Unpacks a single CIGAR string and generates a mapping dictionary.
        Each position is mapped to its corresponding reference coordinate and
        returned before querying"""

        cigar_operations: list = parse_cigar(cigar_string)
        self.check_for_unsupported_ops(cigar_operations)

        mapped_positions = {}
        ref_coord = -1+ref_position
        aln_coord = -1

        for op in cigar_operations:
            for n in range(op.size):

                if op.type in ['M', 'X']:
                    ref_coord, aln_coord = map(lambda x: x + 1, [ref_coord, aln_coord])
                    mapped_positions[aln_coord] = ref_coord
                elif op.type == 'I':
                    aln_coord += 1
                    mapped_positions[aln_coord] = ref_coord
                elif op.type == 'D':
                    ref_coord += 1
                else:
                    raise Exception('Something went wrong during coordinate mapping.  Please check input.')

        return mapped_positions

    def get_all_positions(self, df: pd.DataFrame) -> dict:
        """Fetches all position mappings for a mapping input file and places them into new col in df."""
        all_positions = {}
        for index, row in df.iterrows():
            map: dict = self._map_positions(row['cigar_string'], row['ref_pos'])
            all_positions[row['transcript_id']] = map


        return all_positions


class QueryHandler:
    def __init__(self, map_file, query_file):

        self.input = self._load_map_file(map_file)
        self.query = self._load_query_file(query_file)
        self.mapped_positions = self._generate_maps()

    def _load_map_file(self, file: str) -> pd.DataFrame:
        with open(file, 'r') as f:
            map_df = pd.read_csv(f, sep='\t', names=df_colnames)
            return map_df

    def _load_query_file(self, file: str) -> list:
        with open(file, 'r') as f:
            query_list = []
            for line in f:
                row = line.rstrip('\n').split('\t')
                query_list.append((row[0], int(row[1])))
            return query_list

    def _generate_maps(self) -> dict:
        mapped = CigarMapper(self.input)
        return mapped.position_mappings


def format_output(df: pd.DataFrame):
    """Formats output and saves to a tab-separated file"""
    outfile = df.to_csv(out_file_name, index=False, header=False, sep='\t')
    print(f'Results file: {out_file_name}')

def main(arguments: argparse.Namespace) -> None:

    cigar_in: str = arguments.cigar_path
    query_in: str = arguments.query_path

    query_obj = QueryHandler(cigar_in, query_in)
    input_data: pd.DataFrame =query_obj.input.set_index('transcript_id')
    mappings: dict = query_obj.mapped_positions
    queries: list = query_obj.query

    output = pd.DataFrame(columns=['transcript', 'start', 'chr', 'out'])
    for query in queries:
        transcript, position = query
        if transcript not in mappings.keys():
            line = pd.Series({'transcript': transcript, 'start': position, 'chr':'na', 'out':'transcript not found'})
            output = output.append(line, ignore_index=True)
        else:
            position_map: dict = mappings[transcript]
            chromosome = input_data.at[transcript, 'chromosome']
            try:
                out_coordinate: int = position_map[position]
                line = pd.Series({'transcript': transcript, 'start': position, 'chr': chromosome, 'out': out_coordinate})
                output = output.append(line, ignore_index=True)
            except KeyError:
                line = pd.Series({'transcript': transcript, 'start': position, 'chr': chromosome, 'out': 'position out of range'})
                output = output.append(line, ignore_index=True)

    format_output(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='interface for submitting input files')
    parser.add_argument('--cigar_path', type=str)
    parser.add_argument('--query_path', type=str)
    args = parser.parse_args()

    main(args)



import pytest
import pandas as pd
from cigar_mapper import UnsupportedOperation, CigarMapper, QueryHandler, df_colnames

correct_mappings = {'TR1': {0: 3, 1: 4, 2: 5, 3: 6, 4: 7, 5: 8, 6: 9, 7: 10, 8: 18, 9: 19, 10: 20, 11: 21, 12: 22, 13: 23, 14: 23, 15: 23, 16: 24, 17: 25, 18: 37, 19: 38, 20: 39, 21: 40, 22: 41, 23: 42, 24: 43},
                    'TR2': {0: 10, 1: 11, 2: 12, 3: 13, 4: 14, 5: 15, 6: 16, 7: 17, 8: 18, 9: 19, 10: 20, 11: 21, 12: 22, 13: 23, 14: 24, 15: 25, 16: 26, 17: 27, 18: 28, 19: 29},
                    'TR3': {0: 7, 1: 8, 2: 9, 3: 10, 4: 11, 5: 12, 6: 13, 7: 14, 8: 15, 9: 16, 10: 22, 11: 23, 12: 24, 13: 25, 14: 26, 15: 27, 16: 28, 17: 29, 18: 33}}

@pytest.fixture
def setup_query_obj():
    def _infiles(cigar, query):
        handler = QueryHandler(cigar, query)
        return handler
    return _infiles

@pytest.fixture
def setup_cigar_df():
    def _infile(file):
        with open(file, 'r') as f:
            df =pd.read_csv(file, sep='\t', names=df_colnames)
        return df
    return _infile

@pytest.fixture
def sample_cigar_df(setup_cigar_df):
    df = setup_cigar_df('test_data/test_cigar_in.txt')
    return df


def test_correct_cigars(sample_cigar_df):
    mapper = CigarMapper(sample_cigar_df)
    assert type(mapper) is CigarMapper

def test_invalid_cigar_input(setup_cigar_df):
    df = setup_cigar_df('test_data/test_invalid_cigar.txt')
    with pytest.raises(ValueError):
        mapper = CigarMapper(df)

def test_unsupported_cigar_input(setup_cigar_df):
    df = setup_cigar_df('test_data/test_unsupported_cigar.txt')
    with pytest.raises(UnsupportedOperation):
        mapper = CigarMapper(df)

def test_mappings_match_expected(sample_cigar_df):
    mapper = CigarMapper(sample_cigar_df)
    position_maps: dict = mapper.position_mappings
    for id, map in correct_mappings.items():
        test_map = position_maps[id]
        assert test_map == map

def test_query_handler_setup(setup_query_obj):
    test_handler = setup_query_obj('test_data/test_cigar_in.txt', 'test_data/test_query_in.txt')
    assert type(test_handler) is QueryHandler
    assert type(test_handler.input) is pd.DataFrame
    assert type(test_handler.query) is list
    assert type(test_handler.mapped_positions) is dict
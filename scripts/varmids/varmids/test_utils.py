import pytest
import numpy as np

from varmids.utils import *



# @pytest.fixture
# def known_skew_contents():
#     return {
        
#     }

# def test_nuc_count_calculator(known_skew_contents):
#     for params, vals in known_skew_content.items():
#         nuc_combo = nuc_count_calculator(*params)
#         assert nuc_combo == vals


@pytest.fixture
def known_skews():
    return {
        (20, 20): 0,
        (200, 200): 0,
        (150, 250): 0.25,
        (250, 150): -0.25,
        (1500, 2500): 0.25,
        (2500, 1500): -0.25,
        (0, 100): 1,
        (100, 0): -1
    }


@pytest.fixture
def known_contents():
    kc = {}
    for _ in range(0, 100):
        count_a, count_b = np.random.randint(0, 1000, 2)
        length = np.random.randint(count_a + count_b, count_a + count_b + 1000)
        content = (count_a + count_b) / length
        kc[((count_a, count_b), length)] = content
    return kc


@pytest.fixture
def known_occupied():
    ko = {}
    
    for _ in range(0, 100):
        range_len = np.random.randint(1, 11)
        seq_len = np.random.randint(15, 500)
        sites = np.random.choice([0, 1], seq_len, p=[0.75, 0.25])
        s = np.random.randint(0, seq_len-range_len)
        e = s + range_len
        if sum(sites[s:e]) == 0:
            answer = True
        else:
            answer = False
        ko[answer] = (sites, s, e)
    
    return ko

    
@pytest.fixture
def known_range_lengths():
    rl = {}
    for _ in range(0, 100):
        seq_len = np.random.randint(101, 1000)
        range_len = np.random.randint(1, 100)
        rl[(seq_len, range_len)] = range_len
    return rl


def test_calculate_skew(known_skews):
    for params, val in known_skews.items():
        calculated_skew = calculate_skew(params)
        assert calculated_skew == val, f'{params}: {val} != {calculated_skew}'


def test_calculate_content(known_contents):
    for params, val in known_contents.items():
        calculated_content = calculate_content(*params)
        assert calculated_content == val


def test_random_range_of_length_n(known_range_lengths):
    for params, val in known_range_lengths.items():
        s, e = random_range_of_length_n(*params)
        assert e <= params[0]
        assert e > s
        assert e - s == val


def test_range_is_occupied(known_occupied):
    for val, params in known_occupied.items():
        is_occupied = range_is_occupied(*params)
        assert is_occupied == val







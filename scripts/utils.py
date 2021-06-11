import numpy as np


def range_is_occupied(occupied_coords, start, end):
    sites = np.arange(start, end)
    if any(np.take(occupied_coords, sites)) != 0:
        return True  # is occupied
    else:
        return False


def random_range_of_length_n(length_seq, range_length):
    start = int(np.random.choice(np.arange(0, length_seq-range_length), 1)[0])
    end = start + range_length
    return start, end


def find_available_random_range(occupied_coords, range_length):
    if longest_unoccupied_gap(occupied_coords) >= range_length:
        while True:
            start, end = random_range_of_length_n(
                len(occupied_coords), range_length
                )
            if range_is_occupied(occupied_coords, start, end):
                continue
            else:
                return start, end
    else:
        # no possible ranges
        return False


def longest_unoccupied_gap(occupied_coords):
    # helper function to find longest gap (run of false values) in a boolean
    # array. Used to determine is there is still space in a list for another
    # cluster
    return len(max(''.join([str(i) for i in occupied_coords]).split('1'), key=lambda s: len(s)))
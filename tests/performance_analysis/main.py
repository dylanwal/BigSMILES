from __future__ import annotations

import sys
import time
import datetime
import platform
from collections import deque


import bigsmiles

ZERO_DEPTH_BASES = (str, bytes, int, float, bytearray)


def total_size(object_, verbose: bool = False) -> int:
    """
    Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

    Parameters
    ----------
    object_: Any
        Object you want to get the size of
    verbose: bool
        print out results

    Returns
    -------
    size: int
        size of object in bytes

    """
    seen = set()                      # track which object id's have already been seen
    default_size = sys.getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(obj, name: str):
        if id(obj) in seen:  # do not double count the same object
            return 0

        seen.add(id(obj))
        size = sys.getsizeof(obj, default_size)

        if verbose:
            print(size, name)

        if isinstance(obj, ZERO_DEPTH_BASES):
            pass  # bypass remaining control flow and return
        elif isinstance(obj, (tuple, list, set, deque)):
            size += sum(sizeof(item, name + f"_{i}") for i, item in enumerate(obj))
        elif isinstance(obj, dict):
            size += sum(sizeof(k, name + f".{k}") + sizeof(v, name + f".key{i}")
                        for i, (k, v) in enumerate(getattr(obj, 'items')()))

            # Check for custom object instances - may subclass above too
        if hasattr(obj, '__dict__'):
            size += sizeof(vars(obj), "  " + name + ".")
        if hasattr(obj, '__slots__'):  # can have __slots__ with __dict__
            size += sum(sizeof(getattr(obj, s), "  " + name + f".{s}") for s in obj.__slots__ if hasattr(obj, s))

        return size

    return sizeof(object_, name=type(object_).__name__)


def time_bigsmiles_parsing(polymers: list[str], iter_: int) -> float:
    print(f"Starting time test: {datetime.datetime.now().time()}  (May take a minute.)")
    start_time = time.perf_counter()

    for i in range(iter_):
        bigsmiles.BigSMILES(polymers[i % len(polymers)])

    run_time = time.perf_counter() - start_time
    print(f"Done time test:{datetime.datetime.now().time()}")
    return run_time/iter_ * 1000  # micro-seconds


def time_bigsmiles_parsing_graph(polymers: list[str], iter_: int) -> float:
    print(f"Starting graph time test: {datetime.datetime.now().time()}  (May take a minute.)")
    start_time = time.perf_counter()

    for i in range(iter_):
        bigsmiles.BigSMILES(polymers[i % len(polymers)]).graph()

    run_time = time.perf_counter() - start_time
    print(f"Done graph time test:{datetime.datetime.now().time()}")
    return run_time/iter_ * 1000  # micro-seconds


def memory_bigsmiles_parsing(polymers: list[str]) -> int:
    print(f"Starting memory test: {datetime.datetime.now().time()}")
    data = [bigsmiles.BigSMILES(polymer) for polymer in polymers]
    print(f"Done memory test:{datetime.datetime.now().time()}")
    return int(total_size(data)/len(polymers))  # bytes


def memory_bigsmiles_graph(polymers: list[str]) -> int:
    print(f"Starting memory test: {datetime.datetime.now().time()}")
    data = [bigsmiles.BigSMILES(polymer).graph() for polymer in polymers]
    print(f"Done memory test:{datetime.datetime.now().time()}")
    return int(total_size(data)/len(polymers))  # bytes


def print_memory_breakdown():
    polymer = "[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]"
    result = bigsmiles.BigSMILES(polymer)
    total_size(result, verbose=True)


def main():
    polymer_string = [
        "CC{[>][<]CC(C)[>][<]}CC(C)=C",
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        "CC{[>][$]CC[$],[$]CC(CC)[$][<]}O",
        "CC{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}O",
        "CC{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}F",
        "C{[$][$]CC[$],[$]CC(CC)[$][$]}CC",
        "O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}CC"
    ]
    time_iter = 20_000
    memory_iter = 1000

    bigsmiles_parse_time = time_bigsmiles_parsing(polymer_string, time_iter)
    bigsmiles_parse_memory = memory_bigsmiles_parsing(polymer_string)

    bigsmiles_graph_time = 0 #time_bigsmiles_parsing_graph(polymer_string, time_iter)
    bigsmiles_graph_memory = 0# memory_bigsmiles_graph(polymer_string)

    python_ = sys.version_info
    titles = "date/time (UTF), package version, time per parse (us), memory usage per bigsmiles (bytes), " \
             "platform.processor, graph time (us), memory usage per graph (bytes)", "python version", "notes"
    row = f"\n{datetime.datetime.utcnow()}, {bigsmiles.__version__}, {bigsmiles_parse_time:2.5f}, " \
          f"{bigsmiles_parse_memory:2.0f}, {platform.processor().replace(',', '')}, {bigsmiles_graph_time:2.5f}, " \
          f"{bigsmiles_graph_memory:2.0f}, {python_.major}.{python_.minor}.{python_.micro}, "

    print(titles, row)
    with open("performance.csv", "a", encoding="UTF-8") as file:
        # file.write(titles)
        file.write(row)


if __name__ == "__main__":
    main()
    # print_memory_breakdown()

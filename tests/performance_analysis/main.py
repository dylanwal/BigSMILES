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

    def sizeof(obj):
        if id(obj) in seen:  # do not double count the same object
            return 0

        seen.add(id(obj))
        size = sys.getsizeof(obj, default_size)

        if verbose:
            print(size, type(obj), repr(obj))

        if isinstance(obj, ZERO_DEPTH_BASES):
            pass  # bypass remaining control flow and return
        elif isinstance(obj, (tuple, list, set, deque)):
            size += sum(sizeof(i) for i in obj)
        elif isinstance(obj, dict):
            size += sum(sizeof(k) + sizeof(v) for k, v in getattr(obj, 'items')())

            # Check for custom object instances - may subclass above too
        if hasattr(obj, '__dict__'):
            size += sizeof(vars(obj))
        if hasattr(obj, '__slots__'):  # can have __slots__ with __dict__
            size += sum(sizeof(getattr(obj, s)) for s in obj.__slots__ if hasattr(obj, s))

        return size

    return sizeof(object_)


def time_bigsmiles_parsing(polymers: list[str], iter_: int) -> float:
    print(f"Starting time test: {datetime.datetime.now().time()}  (May take a minute.)")
    start_time = time.perf_counter()

    for i in range(iter_):
        bigsmiles.BigSMILES(polymers[i % len(polymers)])

    run_time = time.perf_counter() - start_time
    print(f"Done time test:{datetime.datetime.now().time()}")
    return run_time/iter_ * 1000  # micro-seconds


def memory_bigsmiles_parsing(polymers: list[str], iter_: int) -> int:
    print(f"Starting memory test: {datetime.datetime.now().time()}")

    data = []
    for i in range(iter_):
        data.append(bigsmiles.BigSMILES(polymers[i % len(polymers)]))

    print(f"Done memory test:{datetime.datetime.now().time()}")
    return int(total_size(data) / iter_)  # bytes


def main():
    polymer_string = [
        "CC{[>][<]CC(C)[>][<]}CC(C)=C",
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        "{[>][$]CC[$],[$]CC(CC)[$][<]}",
        "{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}",
        "{[>][<]C(=O)CCCCC(=O)NCCCCCCN[>][<]}",
        "C{[$][$]CC[$],[$]CC(CC)[$][$]}",
        "[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]"
    ]
    time_iter = 30_000
    memory_iter = 1000

    bigsmiles_parse_time = time_bigsmiles_parsing(polymer_string, time_iter)
    bigsmiles_parse_memory = memory_bigsmiles_parsing(polymer_string, memory_iter)

    titles = "date/time (UTF), package version, time per parse (us), memory usage per bigsmiles (bytes), platform.processor"
    row = f"\n{datetime.datetime.utcnow()}, {bigsmiles.__version__}, {bigsmiles_parse_time:2.5f}, " \
          f"{bigsmiles_parse_memory:2.0f}, {platform.processor().replace(',', '')}"

    print(titles, row)
    with open("performance.csv", "a", encoding="UTF-8") as file:
        # file.write(titles)
        file.write(row)


if __name__ == "__main__":
    main()

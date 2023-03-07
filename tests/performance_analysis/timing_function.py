import functools
import time


def timer(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        print("Finished {} in {} secs".format(repr(func.__name__), round(run_time, 5)))
        return value

    return wrapper


import bigsmiles

n = 100000

@timer
def option1(a, n):
    for _ in range(n):
        b = a == a


@timer
def option2(a, n):
    for _ in range(n):
        b = a is a


a = bigsmiles.BigSMILES("CCCC{[$][$]CC[$][$]}N(NSC)C")
option1(a, n)
option2(a, n)
